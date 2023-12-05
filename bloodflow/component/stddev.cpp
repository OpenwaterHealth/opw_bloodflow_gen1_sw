#include <algorithm>
#include <cmath>

#include "stddev.h"

static void AccumulateFrame(const Frame* src, Frame& dst) {
  for (int i = 0, n = dst.width * dst.height; i < n; i++) {
    dst.data[i] += src->data[i];  // ToDo(jfs): Optimize? (screaming for Eigen vectorization)
  }
}

static void ShiftFrame(Frame& dst) {
  for (int i = 0, n = dst.width * dst.height; i < n; i++) {
    dst.data[i] >>= 6;  // i.e., divide by 64
  }
}

// Calculate mean and variance of a frame, optionally skipping bad and/or masked pixels.

static double Mean(const Frame* fr, const std::vector<bool>& mask, int* satp) {
  // use int64_t to not lose precision on large sums
  // int32_t will overflow for 5MP @ 10bit images
  int64_t sum = 0;
  int size = fr->width * fr->height;
  int count = 0, saturated = 0;
  int saturatedVal = (1 << fr->bits) - 1;

  if (mask.empty()) {
    for (int i = 0; i < size; ++i) {
      sum += fr->data[i];
      if (fr->data[i] == saturatedVal) ++saturated;
    }
    count = size;
  } else {
    for (int i = 0; i < size; ++i) {
      if (mask[i]) {
        sum += fr->data[i];
        if (fr->data[i] == saturatedVal) ++saturated;
        ++count;  // ToDo(jfs): This is the same, once calibration is done.
      }
    }
  }

  if (satp) *satp = saturated;
  return (double)sum / count;
}

static double Variance(const Frame* fr, double mean, const std::vector<bool>& mask) {
  double sum = 0;
  int size = fr->width * fr->height;
  if (mask.empty()) {
    for (int i = 0; i < size; ++i) {
      double var = fr->data[i] - mean;
      sum += var * var;
    }
    return sum / size;
  }

  int count = 0;
  for (int i = 0; i < size; ++i) {
    if (mask[i]) {
      double var = fr->data[i] - mean;
      sum += var * var;
      ++count;
    }
  }
  return sum / count;
}

// Show the current mask by blackening masked pixels.

static void ShowMask(Frame* fr, const std::vector<bool>& mask) {
  if (mask.empty()) return;
  for (int i = 0, size = fr->width * fr->height; i < size; ++i) {
    if (!mask[i]) {
      fr->data[i] = 0;
    }
  }
}

inline int sqr(int i) { return i * i; }

void StdDev::SetROI(int roiX, int roiY, int roiR, int width, int height) {
  roi_.x = roiX, roi_.y = roiY, roi_.r = roiR;
  roi_.width = width, roi_.height = height;

  mask_.clear();  // in case we're called from ReinitROI()
  mask_.resize(width * height, true);
  for (int y = 0, i = 0, r2 = sqr(roiR); y < height; y++) {
    for (int x = 0, y2 = sqr(y - roiY); x < width; x++, i++) {
      if (sqr(x - roiX) + y2 > r2) mask_[i] = false;
    }
  }
}

void StdDev::ReinitROI() {
  SetROI(roi_.x, roi_.y, roi_.r, roi_.width, roi_.height);
}

void* StdDev::Exec(void* data) {
  std::lock_guard<std::mutex> lock(mutex_);  // We'll be changing things in here.
  Frame* fr = (Frame*)data;
  Tag* tag = &pool_.Alloc();

  // Keep a running mean of all frames. This is currently only used by the aligner. There's an
  // experimental history of this file which compares frame means to meanMean, to classify frames
  // as light or dark. Incremental formula for updating mean: Mn = Mn-1 + [Xn - Mn-1] / n.
  int saturated = 0;
  double mean = Mean(fr, mask_, &saturated);
  double variance = Variance(fr, mean, mask_);
  tag->rawMean = mean;
  tag->rawVariance = variance;
  meanMean_ = meanMean_ + (mean - meanMean_) / ++frameCount_;
  // Update the means and accumulate frames. When we've got 64 frames, shift and process.
  if (fr->seq % 2 == 0) {
    // Upate rolling average.
    darkMeanAvg_ = darkMeanAvg_ + (mean - darkMeanAvg_) / ++darkFrameCount_;

    // If accumulation is running, accumulate the frame.
    if (darkAccum_.width && darkAccumCount_ < 64) {
      AccumulateFrame(fr, darkAccum_);
      ++darkAccumCount_;

      // If we've got the frames we want, do the processing for calibration.
      if (darkAccumCount_ == 64) {
        ShiftFrame(darkAccum_);

        // Make a list of bad (hot) pixels from the stats of the averaged frame. If there's a
        // mask, punch holes in it. If not, create one.
        double dMean = Mean(&darkAccum_, std::vector<bool>(), nullptr);
        double dVariance = Variance(&darkAccum_, dMean, std::vector<bool>());
        double dStdDev = sqrt(dVariance);
        int size = fr->width * fr->height, bad = 0;
        if (mask_.empty()) {
          mask_.resize(size, true);
        }
        uint16_t threshold = (uint16_t)(dMean + nSigmas_ * dStdDev);
        for (int i = 0; i < size; i++) {
          if (darkAccum_.data[i] >= threshold) {
            mask_[i] = false;
            ++bad;
          }
        }
        nBadPixels_ = bad;
        darkAveragedMean_ = dMean, darkAveragedStdDev_ = dStdDev;
      }
    }
  } else {
    // Update the rolling average.
    ++lightFrameCount_;
    lightMeanAvg_ = lightMeanAvg_ + (mean - lightMeanAvg_) / lightFrameCount_;
    if (lightAccum_.width && lightAccumCount_ < 64) {
      AccumulateFrame(fr, lightAccum_);
      ++lightAccumCount_;
      if (lightAccumCount_ == 64) {
        ShiftFrame(lightAccum_);
        // ToDo(jfs): auto-circle-fitting for ROI goes here
      }
    }
  }

  if (showMask_) {
    ::ShowMask(fr, mask_);
  }

  tag->saturated = saturated;
  fr->AddTag(tag);
  return data;
}

void StdDev::AtExit(void* data) {
  pool_.Free(GetTag((Frame*)data));
}

void StdDev::resize(size_t size) {
  std::lock_guard<std::mutex> lock(mutex_);
  pool_.resize(size);
  ExecNode::resize(size);
}

void StdDev::beginAccumulation(const Frame* config) {
  std::lock_guard<std::mutex> lock(mutex_);
  darkAccumCount_ = lightAccumCount_ = 0;
  darkAccum_ = Frame(config->width, config->height);
  lightAccum_ = Frame(config->width, config->height);
  size_t size = config->width * config->height * sizeof(darkAccum_.data[0]);
  memset(darkAccum_.data, 0, size);
  memset(lightAccum_.data, 0, size);
}

void StdDev::clearAccumulation() {
  std::lock_guard<std::mutex> lock(mutex_);
  darkAccumCount_ = lightAccumCount_ = 0;
  darkAccum_ = Frame();
  lightAccum_ = Frame();
}

bool StdDev::framesAccumulated() {
  std::lock_guard<std::mutex> lock(mutex_);
  // This is called only periodically, so emit an info msg here.
  // ToDO(jfs): Return counts and print in caller.
  printf("framesAccumulated: %d dark, %d light\n", darkAccumCount_, lightAccumCount_);
  return darkAccumCount_ == 64 && lightAccumCount_ == 64;
}
