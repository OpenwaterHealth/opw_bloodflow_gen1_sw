#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "driver.h"
#include "rcam.h"
#include "rcam_fw.h"
#include "time.h"

#ifdef _MSC_VER
#include <windows.h>  // for SetThreadPriority (TODO(jfs): Move this to a Driver virtual.)
#endif

#define COUNT_ONLY_GOOD_FRAMES 1  // TESTING, e.g., w/old camera and weak laptop

#define CHECK_RET(x)                                       \
  {                                                        \
    int ret = x;                                           \
    if (ret != 0) {                                        \
      printf("Gumstick error: %s returned %d\n", #x, ret); \
      return ret;                                          \
    }                                                      \
  }

// limited register map for HM5530/HM5511/HM1063
#define REG(x, y) static const uint16_t(x) = y;
REG(COARSE_INTEGRATION, 0x0202);
REG(LINE_LENGTH_PCK, 0x0342);
REG(COMMAND_UPDATE, 0x0104);
// (DPC_CONTROL, 0x2481);  silence compiler warning
REG(EXPO_SYNC_CONFIG, 0x2220);
REG(MODE_SELECT, 0x0100);
REG(TS_CFG, 0x2230);
REG(TS_OFFSET_SEL, 0x2233);
REG(TS_VALUE, 0x2237);

// swap endianness issues with I2C communication
static inline uint16_t SwapEndian(uint16_t& u) {
  u = (u >> 8) | (u << 8);
  return u;
}

int Rcam::NumCameras() { return Driver::NumDevices(Rcam::PID, true /* recheck */); }

Rcam::~Rcam() {
  Close();
  delete fx3_;
}

int Rcam::Open(int n) {
  CHECK_RET(fx3_->Open(PID, n));

  CHECK_RET(Flash());

  // get parameters from the gumstick
  static const int MAX_RETRIES = 3;
  int retries;
  for (retries = 0; retries < MAX_RETRIES; ++retries) {
    if (fx3_->CmdRead(RQ_PARAM_READ, 0, 0, sizeof(RcamParam),
                      (uint8_t*)(&param_)) == sizeof(RcamParam)) {
      break;
    }
  }
  if (retries == MAX_RETRIES) {
    printf("Could not read parameters from gumstick\n");
    return -1;
  }

  if (strlen(param_.model) == 0 || strncmp(param_.model, "NONE", 20) == 0) {
    printf("No Image Sensor\n");
    return -1;
  }

  // set up the framebuffer
  size_ = param_.width * param_.height;
  fr_cfg_.bits = param_.bits;
  fr_cfg_.width = param_.width;
  fr_cfg_.height = param_.height;
  x0_ = 0;
  x1_ = param_.width;
  y0_ = 0;
  y1_ = param_.height;
  framebuf_.resize(10, fr_cfg_);
  CHECK_RET(fx3_->BulkInBuffers(1024 * 32, 256, 500));

  frame_length_pck_ = ReadCFG(LINE_LENGTH_PCK, 2);
  // this seems to be the only way to ensure the 1063 is in a known state
  WriteCFG(MODE_SELECT, 0x00);
  Component::SleepMs(100);
  WriteCFG(EXPO_SYNC_CONFIG, 0x02);
  WriteCFG(MODE_SELECT, 0x01);
  Component::SleepMs(100);

  user_ = 0;
  frames_ = 0;
  frames_dropped_ = 0;
  err_ = Frame::OKAY;

  SubWindow2Point(0, 0, param_.width, param_.height);
  SetFrameCount(0);
  return 0;
}

void Rcam::Close() {
  Stop();

  // In synchronous mode (IsLeaf() is true) clear out the framebuffer
  // In async mode, wait until all child processes have finished
  if (IsLeaf()) {
    while (framebuf_.PopAvailable()) framebuf_.Pop();
  } else {
    while (!IsExecDone()) {
      std::this_thread::yield();
    }
  }
}

// The CX3 has difficulty with line widths that are not a multiple of 12.
// The bit packing algorithm causes random USB disconnects, and should not be
// used. Given that the Himax chips do not speed up readout based on x-size, we
// will receive all horizontal pixels, and discard unsed pixels in software.
// Decreasing vertical pixels does result in camera speed up.
void Rcam::SubWindow2Point(int x0, int y0, int x1, int y1) {
  assert(x0 >= 0);
  assert(y0 >= 0);
  assert(x1 > x0);
  assert(y1 > y0);
  assert(x1 <= param_.width);
  assert(y1 <= param_.height);
  assert(rx_state_ == RxState::STOPPED);

  WriteCFG(0x0346, y0, 2);       // Y_ADDR_START
  WriteCFG(0x034A, y1 - 1, 2);   // Y_ADDR_END
  WriteCFG(0x0353, y1 - y0, 2);  // Y_OUTPUT_SIZE;

  // Frame length lines is frame lines + 22, from datasheet
  WriteCFG(0x0340, y1 - y0 + 22, 2);  // FRAME_LENGTH_LINES

  WriteCFG(COMMAND_UPDATE, 1);
  WriteCFG(COMMAND_UPDATE, 0);
  Component::SleepMs(1);

  // reset frame buffer for new image size
  fr_cfg_.width = x1 - x0;
  fr_cfg_.height = y1 - y0;
  size_ = fr_cfg_.width * fr_cfg_.height;
  for (Frame& fr : framebuf_.Raw()) {
    fr = fr_cfg_;
  }

  x1_ = x1;
  x0_ = x0;
  y1_ = y1;
  y0_ = y0;
}

int Rcam::Reset() {
  Stop();
  CHECK_RET(fx3_->Reset());
  CHECK_RET(Flash());

  err_ = Frame::OKAY;

  SubWindow2Point(x0_, y0_, x1_, y1_);
  SetFrameCount(frames_);
  return 0;
}

Frame* Rcam::GetConfig() { return &fr_cfg_; }

void Rcam::SetStream(bool stream) {
  WriteCFG(MODE_SELECT, 0x00);
  WriteCFG(EXPO_SYNC_CONFIG, stream ? 0x00 : 0x02);
  WriteCFG(MODE_SELECT, 0x01);
  Component::SleepMs(1);
}

const char* Rcam::Model() { return param_.model; }

void Rcam::resize(size_t n) {
  assert(rx_state_ == RxState::STOPPED);

  framebuf_.resize(n, fr_cfg_);
  ExecNode::resize(n);
}

double Rcam::Temperature() {
  // ToDo(jfs): Mutex around this? Will that cause bad frames?
  if (strncmp(param_.model, "HM5530", RCAM_PARAM_MODEL_TEXT_ZISE) != 0) {
    return nan("");
  }

  // From Himax:
  // 0x2230TS Enable[0]:  TS _enable; [1]:  TS _auto_mode_enable
  // 0x2233TS _offset_sel[0] : 1, from register; [0] = 0, from OTP
  // 0x2237TS _value_readback[7:0]:8b Hex value
  WriteCFG(TS_CFG, 0x3);
  WriteCFG(TS_OFFSET_SEL, 0x01);
  return ReadCFG(TS_VALUE) * 1.5625 - 45.0;
}

void Rcam::SetGain(double gain) {
  uint16_t gain_reg = (uint16_t)round(log2(gain) * 16.0);
  WriteCFG(0x0205, gain_reg);

  WriteCFG(COMMAND_UPDATE, 1);
  WriteCFG(COMMAND_UPDATE, 0);
  Component::SleepMs(1);
}

void Rcam::SetBLC(bool blc_on) {
  WriteCFG(0x2407, blc_on ? 0x43 : 0x46, 1);  // Change to 43 to turn on (default from digital settings file from himax)
  // Note: himax recommends 0x46 to disable BLC
}

int Rcam::SerialNumber() { return fx3_->SerialNumber(); }

int Rcam::GetFrameCount() { return frames_; }

void Rcam::SetFrameCount(int frames) {
  frames_ = frames;
  int ret = fx3_->CmdWrite(RQ_SET_FRAME_COUNT, 0, 0, 4, (uint8_t*)&frames);
  assert(ret == 4);
}

int Rcam::DroppedFrames() { return frames_dropped_; }

int Rcam::GetErrors() {
  int err = err_;
  err_ = Frame::OKAY;
  return err;
}

void Rcam::SetExposure(double exp) {
  uint16_t coarse_integration = (uint16_t)round(exp * param_.px_clk / frame_length_pck_);
  WriteCFG(COARSE_INTEGRATION, coarse_integration, 2);
  WriteCFG(COMMAND_UPDATE, 1);
  WriteCFG(COMMAND_UPDATE, 0);
  Component::SleepMs(1);
}

double Rcam::GetExposure() {
  uint16_t coarse_integration = ReadCFG(COARSE_INTEGRATION, 2);
  return (double)coarse_integration * (double)frame_length_pck_ /
         (double)param_.px_clk;
}

double Rcam::GetConversion() {
  return (double)frame_length_pck_ / double(param_.px_clk);
}

uint32_t Rcam::ReadCFG(uint32_t reg, int bytes) {
  assert(bytes < 3);
  uint16_t retval = 0;
  int result = fx3_->CmdRead(RQ_I2C_READ, CMOS_ADDR, reg, bytes, (uint8_t*)&retval);
  if (result != bytes) {
    printf("Camera %d: CmdRead failure on register 0x%x; returned %d.\n",
        SerialNumber(), reg, result);
  }
  if (bytes == 2) SwapEndian(retval);
  return retval;
}

void Rcam::WriteCFG(uint32_t reg, uint32_t val, int bytes) {
  assert(bytes < 3);
  uint16_t wval = val;
  if (bytes == 2) SwapEndian(wval);
  fx3_->CmdWrite(RQ_I2C_WRITE, CMOS_ADDR, reg, bytes, (uint8_t*)&wval);  // TODO: on failure?
}

int Rcam::Flash() {
  return fx3_->Flash(__fw_len, __fw);
}

int Rcam::Start() {
  if (rx_state_ == RxState::RUNNING) return 0;

  int result = fx3_->BulkInStart();
  if (result != 0) {
    printf("Failed to start Bulk in endpoint: error %d\n", result);
    rx_state_ = RxState::STOPPED;
    return result;
  }

  frames_bad_ = 0;
  temperature_last_ = Temperature();
  rx_state_ = RxState::RUNNING;
  rx_thread_ = std::thread(&Rcam::RxThread, this);
  temperature_thread_ = std::thread(&Rcam::ReadTemperature, this);
  #ifdef _MSC_VER
  SetThreadPriority(rx_thread_.native_handle(), THREAD_PRIORITY_HIGHEST);
  #endif
  return 0;
}

void Rcam::Stop() {
  if (rx_state_ != RxState::STOPPED) {
    rx_state_ = RxState::QUIT;
    rx_thread_.join();
    temperature_thread_.join();
  }
}

void Rcam::ReadTemperature() {
  while (rx_state_ == RxState::RUNNING) {
    temperature_last_ = Temperature();
    Component::SleepMs(temperature_sample_ms_);
  }
}

//For a patch
bool isSerialNumberValid(int sn) {
  return sn > 0 && sn < 10000;
}
#define RCAM_BAD_SERIAL_NUMBER_ATTEMPTS (5)

void Rcam::RxThread() {
  int serial_number = fx3_->SerialNumber();
  //!!! fx3_->SerialNumber() may return trash
  //A patch
#if 1
  if (!isSerialNumberValid(serial_number)) {
    printf("ERR: unexpected serial number");
    int n = RCAM_BAD_SERIAL_NUMBER_ATTEMPTS;
    for (; --n >= 0;) {
      serial_number = fx3_->SerialNumber();
      if (isSerialNumberValid(serial_number)) {
        break;
      }
      else {
        printf("ERR: still unexpected serial number");
      }
    }
    if (n < 0) {
      printf("ERR: unexpected serial number. Cancelling thread");
      return;
    }
  }
#endif
  // send start command
  fx3_->CmdWrite(RQ_MIPI_START, 0, 0);

  int pixels = 0;
  int x = 0;
  Frame* fr = &framebuf_.Next();
  uint16_t* dst = fr->data;
  uint16_t descramble[4];
  bool gaveUSBError = false;
  while (rx_state_ == RxState::RUNNING) {
    int len = 0;
    uint8_t* buf = NULL;

    int ret = fx3_->BulkInData(len, buf);
    if (ret == Driver::USB_DISCONNECT) {
      if (!gaveUSBError) {
        printf("USB disconnect on camera %d\n", serial_number);
        gaveUSBError = true;
      }
      #if 0  // This was for hot-plugging of cameras ...
      continue;
      #else
      fx3_->BulkInStop();
      fx3_->Reattach();
      fx3_->BulkInStart();
      fx3_->CmdWrite(RQ_MIPI_START, 0, 0);
      len = -1; // signal a packet error
      #endif
    }

    if (len <= 0) continue; // no data available

    len -= PACKET_HEADER_LEN;
    // check for packet errors (invalid length or invalid magic number)
    if (len % 5 != 0) len = -1;
    if (buf == NULL) len = -1;
    if (len >= 0 &&
        strncmp((const char*)(buf + PACKET_HEADER_MAGIC_POS),
                PACKET_HEADER_MAGIC,
                PACKET_HEADER_MAGIC_LEN) != 0) {
      len = -1;
    }
    uint8_t* in = buf + PACKET_HEADER_LEN;
    while (len > 0) {
      descramble[0] = ((uint16_t)in[0] << 2) | ((in[4] >> 0) & 3);
      descramble[1] = ((uint16_t)in[1] << 2) | ((in[4] >> 2) & 3);
      descramble[2] = ((uint16_t)in[2] << 2) | ((in[4] >> 4) & 3);
      descramble[3] = ((uint16_t)in[3] << 2) | ((in[4] >> 6) & 3);
      for (int i = 0; i < 4; ++i) {
        if (x >= x0_ && x < x1_) {
          if (dst) *(dst++) = descramble[i];
          ++pixels;
          if (pixels > size_) break;
        }
        ++x;
        if (x >= param_.width) x = 0;
      }
      len -= 5;
      in += 5;
      if (pixels > size_) break;
    }

    // check for EOF
    if ((len < 0) || (buf[PACKET_EOF_POS] & PACKET_EOF_VAL) || (pixels >= size_)) {

      // check for errors
      int err = Frame::OKAY;
      if (len < 0) {
        err |= Frame::ERR_PACKET;
      } else {
        #if COUNT_ONLY_GOOD_FRAMES == 0
        // move device_frames from the packet into device_frames_ atomically
        // buf[FRAME_COUNT_POS] may not be 32-bit aligned, so use memcpy first
        int device_frames_buf;
        memcpy(&device_frames_buf, &buf[FRAME_COUNT_POS], FRAME_COUNT_BYTES);
        if (device_frames_buf != frames_) err |= Frame::ERR_ORDER;
        frames_ = device_frames_buf;
        #endif
      }
      if (pixels > size_) err |= Frame::ERR_OVERFLOW;
      if (pixels < size_) err |= Frame::ERR_INCOMPLETE;
      err_ |= err;

      if (err != Frame::OKAY) { ++frames_bad_; }

      // Only emit a frame if:
      //   1. fr pointer is valid (we were storing data as this frame came in)
      //   2. len == 0 (we had a valid packet header)
      //   3. FRAME_VALID is true
      if (len == 0 && (buf[FRAME_VALID_POS] & FRAME_VALID_VAL)) {
        if (fr) {
          fr->bits = param_.bits;
          fr->width = fr_cfg_.width;
          fr->height = fr_cfg_.height;
          fr->ClearTags();
          fr->serialNumber = serial_number;
          fr->err = err;
          fr->seq = frames_;
          fr->temperature = temperature_last_;
          fr->SetTimestamp();
          framebuf_.Push();
          if (!IsLeaf()) { /*printf("cam %d, frame %d\n", serial_number, fr->seq);*/ Produce(fr); }
          fr = NULL;
          #if COUNT_ONLY_GOOD_FRAMES
          ++frames_;
          #endif
        }
        #if COUNT_ONLY_GOOD_FRAMES == 0
        ++frames_;
        #endif
      }

      if (!fr && framebuf_.PushAvailable()) {
        fr = &framebuf_.Next();
      }

      if (fr) {
        dst = fr->data;
      } else {
        dst = NULL;
        ++frames_dropped_;
      }

      pixels = 0;
      x = 0;
    }
  }

  fx3_->CmdWrite(RQ_MIPI_STOP, 0, 0);
  fx3_->BulkInStop();
  rx_state_ = RxState::STOPPED;
}

void Rcam::AtExit(void* data) {
  #if COUNT_ONLY_GOOD_FRAMES == 0
  if ((Frame*)data != &framebuf_.Peek()) {
    printf("Out of order frame %d, expected %d\n", ((Frame*)data)->seq,
           framebuf_.Peek().seq);
  }
  assert((Frame*)data == &framebuf_.Peek());
  framebuf_.Pop();
  #else
  while (framebuf_.PopAvailable()) framebuf_.Pop();
  #endif
}
