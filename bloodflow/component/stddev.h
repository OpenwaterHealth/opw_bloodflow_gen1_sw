#pragma once

#include <mutex>
#include <vector>

#include "execnode.h"
#include "frame.h"
#include "pool.h"

// Compute the mean and standard deviation on frames
class StdDev : public ExecNode {
 public:
  StdDev() {}
  ~StdDev() {}

  struct Tag : Frame::Tag {
    double rawMean;     // spatial mean of this frame
    double rawVariance; // spatial variance of this frame
    int saturated;      // count of saturated pixels
  };

  // Get a mean / standard deviation tag from a frame
  // returns NULL if no tag found
  // @note deprecated, use Frame::GetTag<StdDev::Tag>();
  static Tag* GetTag(const Frame* frame) { return frame->GetTag<Tag>(); }

  // resize the number of ExecNodes that can be active at one time
  void resize(size_t size) override;

  // Reset dark/light averaging process.
  void resetAveraging() {
    darkFrameCount_ = lightFrameCount_ = 0;
    darkMeanAvg_ = lightMeanAvg_ = 0;
  }

  // Return the running means.
  double meanMean() { return meanMean_; }
  double darkMean() { return darkMeanAvg_; }
  double lightMean() { return lightMeanAvg_; }

  // Frame accumulation averages 64 frames to reduce noise, in memory. We have 10-bit images from
  // the camera, leaving us 6 bits to spare. After which we can divide by 64 easily.

  // Start averaging frames for calibration.
  void beginAccumulation(const Frame* config);

  // Clear the accumulators (notably returning the memory used).
  void clearAccumulation();

  // Return accumulator status.
  bool framesAccumulated();

  // Return accumulated frames.
  Frame& darkAccum() { return darkAccum_; }
  Frame& lightAccum() { return lightAccum_; }

  // Set the beam ROI. Not reentrant; expected to be called at set-up time.
  void SetROI(int x, int y, int r, int width, int height);

  // Reinit the ROI from stored params.
  void ReinitROI();

  // Zero the frame count.
  void ResetFrameCount() { frameCount_ = 0; }

  // Set & get gain constant
  void SetGainConstant(double k) { gainConst_ = k; }
  double GainConstant() const { return gainConst_; }

  // Show the mask (or not).
  void ShowMask(bool showMask) { showMask_ = showMask; }

  // accessors for bad pixel detection params & stats
  double nSigmas() { return nSigmas_; }
  int nBadPixels() { return nBadPixels_; }
  double darkAveragedMean() { return darkAveragedMean_; }
  double darkAveragedStdDev() { return darkAveragedStdDev_; }

 protected:
  void* Exec(void* data) override;
  void AtExit(void* data) override;

  Pool<Tag> pool_;

  // counts of frames contributing to means
  int frameCount_ = 0, darkFrameCount_ = 0, lightFrameCount_ = 0;

  // rolling window of means of all frames, dark frames and light frames
  double meanMean_ = 0, darkMeanAvg_ = 0, lightMeanAvg_ = 0;

  // mean & variance  of most recent dark frame (presumably, immediately before a light frame)
  double darkMean_, darkVar_;

  // Accumulators for averaging, and frame counts thereto.
  Frame darkAccum_, lightAccum_;
  int darkAccumCount_, lightAccumCount_;

  // default # of std-devs above mean to classify bad (hot) pixels
  double nSigmas_ = 3;

  // stats from bad pixel detection
  int nBadPixels_ = 0;
  double darkAveragedMean_ = 0;
  double darkAveragedStdDev_ = 0;

  // circular ROI for beam center
  struct ROI {
    int x = 0, y = 0, r = 0, width = 0, height = 0;
  } roi_;

  double gainConst_ = 0;  // Optional gain constant
  bool showMask_ = false;  // show the mask? (by blackening pixels)

  // ROI as per-pixel mask. We could speed things up a bit by computing an x range for each y, but
  // this is easy, and includes the bad pixel mask. Recall std::vector<bool> is specialized to use
  // only one bit per element, so it uses less memory .. and smaller means faster.
  std::vector<bool> mask_;

  std::mutex mutex_;
};
