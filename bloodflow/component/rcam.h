#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <atomic>

#include "circular_buffer.h"
#include "execnode.h"
#include "frame.h"
#include "rcam_param.h"

class Driver;

// Class to connect to cameras and read raw data
// Contains an internal circular framebuffer to store data as it comes
//   in from the device
// Two APIs exist:
// - Legacy API: Synchonous using GetFrame() or WaitFrame()
//   Rcam considers a frame "read" once GetFrame has been called to read a
//   subsequent frame (even if no subsequent frame is available), and may
//   start writing over the previous frame
// - New API: Asynchronous using ExecNodes to process data
//   Rcam will keep frames as long as subesquent ExecNodes are still
//   operating on the data.
class Rcam : public ExecNode {
 public:
  // check how many Rcams are connected to this computer
  static int NumCameras();

  Rcam(Driver* device) { fx3_ = device; }  // assumes ownership of fx3 (TODO(jfs): smartptr)
  virtual ~Rcam();

  // open a camera
  // @param n which camera, 0 -> N-1 in a system with N cameras
  //   to connect to
  // @returns 0 on success, error code otherwise
  virtual int Open(int n = 0);

  // close the currently open camera
  void Close();

  // Return the camera serial number (via driver).
  virtual int SerialNumber();

  // Start the camera. Returns 0 on success, error code otherwise.
  virtual int Start();

  // Stop the camera
  virtual void Stop();

  // Subwindow the camera, only receiving pixels in the box defined by
  // upper left: (x0, y0) lower right: (x1, y1)
  // x1 and y1 are the first pixels NOT displayed, so
  //   SubWindow2Point(0, 0, width, height) will display the entire frame,
  //   and SubWindow2Point(0, 0, 1, 1) will display nothing
  // @note Call this immediately after Open(), as it will change the
  //   configuration frame.
  virtual void SubWindow2Point(int x0, int y0, int x1, int y1);

  // Hard reset of camera firmware
  // Will call Close() if camera is running, and unload firmware
  // Pauses for 2.5 seconds to allow USB reattach
  // Call Open() to restart communication
  // @returns 0 on success, otherwise an error code
  virtual int Reset();

  // Set the camera to auto-stream at its programmed FPS (ignores SYNC pin)
  virtual void SetStream(bool stream);

  // Get configuration frame
  // Frame info (height, width, bit depth) will be valid
  // Frame specifics (data, timestamp, etc..) will not
  Frame* GetConfig();

  // check number of received frames
  virtual int GetFrameCount();

  // Set the number of received frames.
  // The next frame received will have its sequence number set to this.
  void SetFrameCount(int frames);

  // Check number of dropped frames.
  int DroppedFrames();

  // Check number of bad frames.
  int BadFrames() const { return frames_bad_; }

  // Check for errors since GetErrors() was last called
  // @returns a bitwise OR of all errors that have occurred since last call
  int GetErrors();

  // String describing the camera model (ie, "HM5530")
  const char* Model();

  // Set exposure
  // @param exp exposure time in seconds
  virtual void SetExposure(double exp);

  // Get exposure
  // @returns exposure time in seconds
  double GetExposure();

  // Get conversion time
  // @returns conversion time in seconds
  double GetConversion();

  // Get temperature of CMOS sensor in degrees C
  double Temperature();

  // Set the analog gain
  virtual void SetGain(double gain);

  // Turn on/off black level correction
  virtual void SetBLC(bool blc_on);

  // write configuration to the device
  void WriteCFG(uint32_t reg, uint32_t val, int bytes = 1);

  // read configuration from the device
  uint32_t ReadCFG(uint32_t reg, int bytes = 1);

  // resize the buffer for a number of frames
  void resize(size_t n) override;

  // for testing
  static int pid() { return PID; }

 protected:
  // Underlying FX3 driver
  Driver* fx3_ = nullptr;

  static const int PID = 0x4F10;

  void AtExit(void* data) override;

  void RxThread();
  volatile enum class RxState { STOPPED, RUNNING, QUIT } rx_state_ = RxState::STOPPED;
  std::thread rx_thread_;

  RcamParam param_ = {{0}};
  Frame fr_cfg_;
  CircularBuffer<Frame> framebuf_;

  int Flash();

  // number of frames read, written, and in userspace (0 or 1)
  // in this session
  volatile int user_ = 0;
  volatile int frames_ = 0;
  volatile int frames_dropped_ = 0;
  volatile int frames_bad_ = 0;
  volatile int err_ = Frame::OKAY;

  // conversion time in pixel clocks (exposure time is set in
  // units of conversion times)
  uint16_t frame_length_pck_;

  // track device temperature
  time_t temperature_sample_ms_ = 100;
  std::atomic<double> temperature_last_;
  void ReadTemperature();
  std::thread temperature_thread_;

  // track x subwindow
  int x0_ = 0;
  int x1_ = 0;
  int y0_ = 0;
  int y1_ = 1;
  int size_ = 0;
};
