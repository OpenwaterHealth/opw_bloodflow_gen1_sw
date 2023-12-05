#include <cstring>
#include <assert.h>

#include "bloodflow/third_party/tiff/inc/tiffio.h"
#include "bloodflow/third_party/tiff/inc/tiff.h"

#include "frame.h"
#include "time.h"

// shift right/left based on sign
static inline uint8_t sh(uint16_t in, int shift) {
  if (shift >= 0) {
    return in >> shift;
  } else {
    return in << (-shift);
  }
}

class RealTiff: public TiffInterface {
 public:
  RealTiff() {}
  ~RealTiff() {}

  bool Open(const char* file_name, const char* mode) {
    tiff_ = TIFFOpen(file_name, mode);
    return tiff_ != nullptr;
  }

  void Close() { TIFFClose(tiff_); tiff_ = nullptr; }

  int SetField(ttag_t tag, int value) { return TIFFSetField(tiff_, tag, value); }

  int WriteScanline(tdata_t buf, uint32 row) { return TIFFWriteScanline(tiff_, buf, row); }

  int WriteDirectory() { return TIFFWriteDirectory(tiff_); }

 private:
  TIFF* tiff_;
};

void Frame::Init() {
  data = nullptr;
  bits = 0;
  width = 0;
  height = 0;
  seq = 0;
  serialNumber = -1;
  timestamp_ms_ = 0;

  err = Frame::OKAY;
  tiff_ = nullptr;
}

Frame::Frame() {
  Init();
}

Frame::Frame(int w, int h) {
  Init();
  width = w;
  height = h;
  bits = 16;
  data_.resize(w * h);
  data = data_.data();
}

Frame::Frame(const Frame& fr) {
  bits = fr.bits;
  width = fr.width;
  height = fr.height;
  seq = fr.seq;
  serialNumber = fr.serialNumber;
  timestamp_ms_ = fr.timestamp_ms_;
  tags_ = fr.tags_;
  err = fr.err;
  if (fr.data) {
    data_ = fr.data_;
  } else {
    data_.resize(width * height);
  }
  data = data_.data();
  tiff_ = nullptr;  // Do not copy this!
}

const Frame& Frame::operator=(const Frame& fr) {
  bits = fr.bits;
  width = fr.width;
  height = fr.height;
  seq = fr.seq;
  serialNumber = fr.serialNumber;
  timestamp_ms_ = fr.timestamp_ms_;
  tags_ = fr.tags_;
  err = fr.err;
  if (fr.data) {
    data_ = fr.data_;
  } else {
    data_.resize(width * height);
  }
  data = data_.data();
  tiff_ = nullptr;  // Do not copy this!

  return *this;
}

Frame::~Frame() {
  data = nullptr;
  delete tiff_;
  tiff_ = nullptr;
}

uint16_t* Frame::operator[](int col_idx) {
  return data + col_idx * width;
}

void Frame::SetTimestamp() {
  timestamp_ms_ = Component::SteadyClockTimeMs();
}

void Frame::InitTiff() {
  // If you cannot afford a TiffInterface, one will be appointed for you ...
  if (!tiff_) {
    tiff_ = new RealTiff();
  }
}

int Frame::Write(const std::string& fname) {
  InitTiff();
  // if (!tiff_) return -1;  // never happens
  if (!tiff_->Open(fname.c_str(), "w")) return -1;

  class stackTiff {  // Close on destruction.
   public:
    stackTiff(TiffInterface* tiff): t(tiff) {}
    ~stackTiff() { t->Close(); }
    TiffInterface* t;
  };
  stackTiff st(tiff_);

  tiff_->SetField(TIFFTAG_IMAGEWIDTH, width);
  tiff_->SetField(TIFFTAG_IMAGELENGTH, height);
  tiff_->SetField(TIFFTAG_BITSPERSAMPLE, bits);
  tiff_->SetField(TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  tiff_->SetField(TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  tiff_->SetField(TIFFTAG_ORIENTATION, static_cast<int>(ORIENTATION_TOPLEFT));
  tiff_->SetField(TIFFTAG_SAMPLESPERPIXEL, 1);
  tiff_->SetField(TIFFTAG_ROWSPERSTRIP, 1);
  tiff_->SetField(TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  tiff_->SetField(TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
  tiff_->SetField(TIFFTAG_MINSAMPLEVALUE, 0);
  tiff_->SetField(TIFFTAG_MAXSAMPLEVALUE, (1 << bits) - 1);

  if (bits == 16) {
    // For Rcam, this will write 10-bit data as 16, just wasting some space.
    for (int j = 0; j < height; ++j) {
      if (tiff_->WriteScanline((*this)[j], j) != 1) return -1;
    }
  } else if (bits == 8) {
    // For assumed 10-bit data, we want the top 8 bits, not the bottom.
    uint8_t* line_data = new uint8_t[width];
    for (int j = 0; j < height; ++j) {
      for (int i = 0; i < width; ++i) {
        line_data[i] = sh(data[i + width * j], 2);
      }
      if (tiff_->WriteScanline(line_data, j) != 1) return -1;
    }
    delete[] line_data;
  } else {
    int line_len = width * bits / 8;
    if (width * bits % 8 != 0) ++line_len;
    uint8_t* line_data = new uint8_t[line_len];
    for (int j = 0; j < height; ++j) {
      uint8_t* line_ptr = line_data;
      memset(line_data, 0, line_len);
      int shift = 0;
      for (int i = 0; i < width; ++i) {
        shift += bits;
        while (true) {
          *line_ptr |= sh(data[i + width * j], shift - 8);
          if (shift >= 8) {
            ++line_ptr;
            shift -= 8;
            if (shift == 0) break;
          } else {
            break;
          }
        }
      }
      if (tiff_->WriteScanline(line_data, j) != 1) return -1;
    }
    delete[] line_data;
  }

  if (tiff_->WriteDirectory() != 1) return -1;
  return 0;  // stackTiff closes
}
