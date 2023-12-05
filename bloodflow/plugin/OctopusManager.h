#pragma once

#include "bloodflow/component/octopus.h"
#include "bloodflow/third_party/json-develop/single_include/nlohmann/json.hpp"

// Encapsulate all the expertise on octopus timers for scanning
class OctopusManager {
 public:
  OctopusManager(Octopus* octopus): octopus_(octopus) {}  // takes ownership of octopus
  virtual ~OctopusManager();

  using json = nlohmann::json;

  // Establish usb connection to octopus, calls InitializeOctopus function
  virtual bool init(const json& systemParameters);

  // Software trigger to initiate software gate to Frame Valid and AOM Chopper
  virtual void TriggerDataCollection();

  // Enable output from software and hardware triggers set up in InitializeOctopus function
  virtual bool EnableSystemChannels(bool channelON);

  // Configure channels after reset, disable, or end state reached. This will enable multiple data
  // acquisitions of a fixed number of frames (transitions).
  virtual bool ConfigureSystemChannels();

  // Enable frame valid output [set channel to high/low to ensure high before FSIN signals start]
  virtual bool EnableFrameValid(bool validHIGH);

  // Return the Octopus.
  Octopus* GetOctopus() { return octopus_; }

 protected:
  // Set up software trigger/gate/clock, and hardware triggers to FSIN, frame valid, and AOM chopper
  bool InitializeOctopus(const json& systemParameters);

  // Set pin used to trigger data acquisition (hardware vs software)
  bool SetDataCollectionTriggerPin(Octopus::Pin triggerPin);

  // Set up AOM chopper for pseudo-pulsed systems; used to change pulsewidths during scanning
  // (currently used only internally, but called from outside for multiple exposure system)
  bool SetAOMPulsewidth(int numImages, double delay, double width, double freq_Hz, double volt_V);

  // Set all system timers to low
  bool DisableOutput(Octopus::Timer timer, Octopus::Pin pin);

  Octopus* octopus_;
  struct systemTimer {
    Octopus::Timer timer;
    Octopus::Pin pin;
  };

  // Timers/pins that are for internally triggering and gating other stuff
  Octopus::Timer timer_SWtrigger_;  // SW trigger to initiate data collection collection [for systems with no pushbutton switch]
  systemTimer swClock_;
  systemTimer swGate_;
  systemTimer wandLED_;
  systemTimer aomPulseChopGate_;  // SW gate for pulse chopping AOM time
  Octopus::Pin pin_SWgate_ = Octopus::Pin::OUT8_TOP;  // SW gate to frame valid and AOM chopper [if hw triggered system this is a hw gate]
  Octopus::Pin pin_SWclock_ = Octopus::Pin::OUT7_TOP; // SW clock to sync camera and AOM chopper [option for future to replace with hw input to sync with laser]
  Octopus::Pin pin_AOMpulseChopGate_ = Octopus::Pin::OUT5_TOP;
  int channel_AOMpulseChop_ = 4;  // analog output channel for AOM pulse chopper
  int hwTriggered_ = 0;  // is system hw/button or sw triggered?

  // Camera Pins FSIN & Frame Valid trigger pins
  std::vector<Octopus::Pin> pins_cameraFSIN_ = {
      Octopus::Pin::OUT1_BOTTOM, Octopus::Pin::OUT2_BOTTOM,
      Octopus::Pin::OUT3_BOTTOM, Octopus::Pin::OUT4_BOTTOM };
  std::vector<Octopus::Pin> pins_cameraFrameValid_ = {
      Octopus::Pin::OUT1_TOP, Octopus::Pin::OUT2_TOP,
      Octopus::Pin::OUT3_TOP, Octopus::Pin::OUT4_TOP };
  systemTimer tempFrameValidTimer_;
  std::vector<systemTimer> timers_cameraFrameValid_;
  std::vector<systemTimer> timers_cameraFSIN_;

  // LED I/O
  Octopus::Pin pin_wandLED_ = Octopus::Pin::OUT6_BOTTOM; // Illuminates when "okay" for data collection, shuts off when data collection finished

  // Data collection and syncing related pins
  // Trigger for data collection
  // Can be hardware (push button switch TTL input) or software. Configued using SetDataCollectionTriggerPin
  Octopus::Pin pin_dataCollectionTrigger_ = Octopus::Pin::LOW; // Default to no triggers for safety
};
