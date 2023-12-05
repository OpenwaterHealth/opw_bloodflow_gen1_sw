#pragma once

#include <mutex>

// CyAPI has issues if windows.h is not included first
#include <windows.h>
#include "bloodflow/third_party/cypress-fx3/inc/CyAPI.h"

#include "driver.h"

// Interface to generic FX3 devices.
class FX3: public Driver {
 public:
  FX3() {}
  ~FX3() { Close(); }

  // Connect to a device with a particular PID and serial number
  // @param pid PID of the device
  // @param n nth device with the specified PID to which to connect
  //   Leave blank to connect to any device found.
  // @returns 0 on success
  int Open(uint16_t pid, int n = 0) override;

  // Reset the target device to bootloader
  int Reset() override;

  // Serial number of this device
  // @returns serial number if found, -1 otherwise
  int SerialNumber() override;

  // Close connection to the FX3 device.
  void Close() override;

  // Flash an FX3 with a firmware image
  // @param len length of the firmware image
  // @param data data to flash
  // @returns 0 on success
  int Flash(int len, uint8_t* data) override;

  // Threadsafe write control packet to endpoint 0
  // @param req request code
  // @param wValue request value
  // @param wIndex request index
  // @param len length of request, if blank 0
  // @param buf bytes with length len to transmit
  // @returns bytes written
  int CmdWrite(uint8_t req, uint16_t wValue, uint16_t wIndex, uint16_t len = 0, uint8_t* buf = NULL)
      override;

  // Threadsafe read control packet to endpoint 0
  // @param req request code
  // @param wValue request value
  // @param wIndex request index
  // @param len length of request
  // @param buf bytes with length len to fill with recieved data
  // @returns bytes read
  int CmdRead(uint8_t req, uint16_t wValue, uint16_t wIndex, uint16_t len, uint8_t* buf) override;

  // Receive data from a bulk in endpoint in synchronous mode
  // @param len maximum number of bytes to receive
  // @param data buffer to fill with data
  // @returns number of bytes read, <0 if an error has occurred
  int DataIn(int len, uint8_t* data) override;

  // Transmit data to a bulk in endpoint in synchronous mode
  // @param len number of bytes to transmit
  // @param data data to transmit
  // @returns number of bytes transmitted, <0 if an error has occurred.
  int DataOut(int len, uint8_t* data) override;

  // Bulk in asynchronous control
  // The bulk in endpoint operates in blocking mode by default. This is sufficient for low data
  // rate endpoints, or synchronous endpoints. For high data rate endpoints, it is recommended to
  // use the BulkInBuffers() and BulkInStart() functions to setup the buffers and start streaming
  // data in respectively. BulkInData() will then draw data from those buffers as they fill.
  // BulkInStop() is automatically called by the destructor, but can be explicitly called as well.

  // Set buffers for bulk in endpoint
  // @param buflen buffer length for a single bulk transaction
  //   if 0, will revert to synchronous mode
  // @param n_buffers number of buffers per endpoint.
  // @param timeout_ms milliseconds to wait for each buffer to be filled
  // @returns 0 on success
  int BulkInBuffers(int buflen, int n_buffers, int timeout_ms) override;

  // Start bulk input transfer
  // @returns 0 on success
  int BulkInStart() override;

  // Get data from the bulk in endpoint
  // @param[out] len length of data produced
  // @param[out] data pointer to data
  // @returns 0 on success
  // @note may block if waiting on data
  int BulkInData(int& len, uint8_t*& data) override;

  // Stop bulk input transfer
  // @returns 0 on success
  void BulkInStop() override;

  // Reattach the device and refresh endpoint pointers
  void Reattach() override;

 protected:
  static const int RQ_SERIAL = 0xE1;
  static const int RQ_RESET = 0xE0;

  // device endpoints
  CCyUSBDevice usb_device_;
  CCyControlEndPoint* ctrl_ = NULL;
  CCyBulkEndPoint* bulk_in_ = NULL;
  CCyBulkEndPoint* bulk_out_ = NULL;
  int address_ = -1;

  // bulk in buffers
  bool async_ = false;
  int timeout_ms_ = 0;
  std::vector<std::vector<uint8_t>> in_buffer_;
  std::vector<OVERLAPPED> inovlap_;
  std::vector<uint8_t*> contexts_;
  int in_n_ = -1;

  std::mutex ctrl_mutex_;
};
