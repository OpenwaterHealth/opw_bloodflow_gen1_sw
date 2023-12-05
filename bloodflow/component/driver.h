#pragma once

#include <cstdint>
#include <vector>

// Pure-virtual Base class for interface to generic FX3 devices.
class Driver {
 public:
  Driver() {}
  virtual ~Driver() {}

  // Error codes
  // USB disconnect error
  static const int USB_DISCONNECT = -1000;

  // Number of FX3 with a specific PID attached to the computer.
  // @note PID for Openwater devices have the lower bit cleared
  //   when unprogrammed, and set when programmed.
  //   This means PID 0x1A0 and 0x1A1 are the same physical device,
  //   however, PID 0x1A1 is programmed and ready for instructions.
  //   This was chosen to allow programmed / unprogrammed devices
  //   to be visible from the device manager.
  // N.B.: Derived classes must define this as Driver::, since it's static.
  static int NumDevices(uint16_t pid, bool recheck = false);

  // Connect to a device with a particular PID and serial number
  // @param pid PID of the device
  // @param n nth device with the specified PID to which to connect
  //   Leave blank to connect to any device found.
  // @returns 0 on success
  virtual int Open(uint16_t pid, int n = 0) = 0;

  // Reset the target device to bootloader
  virtual int Reset() = 0;

  // Serial number of this device
  // @returns serial number if found, -1 otherwise
  virtual int SerialNumber() = 0;

  // Close connection to the FX3 device.
  virtual void Close() = 0;

  // Flash an FX3 with a firmware image
  // @param len length of the firmware image
  // @param data data to flash
  // @returns 0 on success
  virtual int Flash(int len, uint8_t* data) = 0;

  // Threadsafe write control packet to endpoint 0
  // @param req request code
  // @param wValue request value
  // @param wIndex request index
  // @param len length of request, if blank 0
  // @param buf bytes with length len to transmit
  // @returns bytes written
  virtual int CmdWrite(uint8_t req, uint16_t wValue, uint16_t wIndex, uint16_t len = 0,
                       uint8_t* buf = nullptr) = 0;

  // Threadsafe read control packet to endpoint 0
  // @param req request code
  // @param wValue request value
  // @param wIndex request index
  // @param len length of request
  // @param buf bytes with length len to fill with recieved data
  // @returns bytes read
  virtual int CmdRead(uint8_t req, uint16_t wValue, uint16_t wIndex, uint16_t len,
                      uint8_t* buf) = 0;

  // Receive data from a bulk in endpoint in synchronous mode
  // @param len maximum number of bytes to receive
  // @param data buffer to fill with data
  // @returns number of bytes read, <0 if an error has occurred
  virtual int DataIn(int len, uint8_t* data) = 0;

  // Transmit data to a bulk in endpoint in synchronous mode
  // @param len number of bytes to transmit
  // @param data data to transmit
  // @returns number of bytes transmitted, <0 if an error has occurred.
  virtual int DataOut(int len, uint8_t* data) = 0;

  // Bulk in asynchronous control
  // The bulk in endpoint operates in blocking mode by default.
  // This is sufficient for low data rate endpoints, or synchronous endpoins.
  // For high data rate endpoints, it is recommended to use the
  //   BulkInBuffers() and BulkInStart() functions to setup
  //   the buffers and start streaming data in respectively.
  //   BulkInData() will then draw data from those buffers as they
  //   fill.  BulkInStop() is automatically called by the destructor, but can
  //   be explicitly called as well.

  // Set buffers for bulk in endpoint
  // @param buflen buffer length for a single bulk transaction
  //   if 0, will revert to synchronous mode
  // @param n_buffers number of buffers per endpoint.
  // @param timeout_ms milliseconds to wait for each buffer to be filled
  // @returns 0 on success
  virtual int BulkInBuffers(int buflen, int n_buffers, int timeout_ms) = 0;

  // Start bulk input transfer
  // @returns 0 on success
  virtual int BulkInStart() = 0;

  // Get data from the bulk in endpoint
  // @param[out] len length of data produced
  // @param[out] data pointer to data
  // @returns 0 on success
  // @note may block if waiting on data
  virtual int BulkInData(int& len, uint8_t*& data) = 0;

  // Stop bulk input transfer
  // @returns 0 on success
  virtual void BulkInStop() = 0;

  // Reattach the device and refresh endpoint pointers
  virtual void Reattach() = 0;

 protected:
  // These are defined here for subclasses to define. See fx3.cpp for the real versions. Test
  // modules can skip these, needing only NumDevices().

  // The Cypress driver maps Open(device_number) to different devices after a USB
  //   disconnect.  To work around this, we use the underlying windows ID, and map
  //   that to device number.
  // Because of this, we need to create a list of device PIDs and addresses the first
  //   time this driver is called, to ensure Open(0) always opens the 0th device
  //   connected to the system (ie, the cypress driver hasn't reordered devices between
  //   Open(0); ... Open(1);
  static void Init(bool recheck = false);

  // map of PID to address.  may have multiple idendical entries for PID, which
  //   correspond to multiple devices of the same type connected to the system.
  static std::vector<std::pair<uint16_t, int>> addresses_;

  // Find the device number that has address.  Device numbers change based on USB
  //   disconnects, addresses do not.
  static int DeviceNumber(int address);
};
