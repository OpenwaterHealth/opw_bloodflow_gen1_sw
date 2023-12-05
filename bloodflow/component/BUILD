package(default_visibility = ["//visibility:public"])

config_setting(name = "win", constraint_values = ["@platforms//os:windows"])
config_setting(name = "mac", constraint_values = ["@platforms//os:osx"])
config_setting(name = "linux", constraint_values = ["@platforms//os:linux"])

lcov_copts = select({
  ":linux": [ "-g -O0 -fprofile-instr-generate -fcoverage-mapping" ],
  "//conditions:default": [],
})
lcov_linkopts = select({
  ":linux": [ "-g -O0 -fprofile-instr-generate -fcoverage-mapping" ],
  "//conditions:default": [],
})

cc_library(
  name = "circular_buffer",
  hdrs = [ "circular_buffer.h" ],
)

cc_library(
  name = "driver",
  hdrs = [ "driver.h" ],
)

cc_library(
  name = "execnode",
  hdrs = [ "execnode.h" ],
  srcs = [ "execnode.cpp" ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [
    ":circular_buffer",
  ],
)

cc_test(
  name = "execnode_test",
  srcs = [ "test/execnode_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":execnode",
    ":time",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ],
)

cc_library(
  name = "frame",
  hdrs = [ "frame.h" ],
  srcs = [ "frame.cpp" ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [
    ":tiff_interface",
    ":time",
    "//bloodflow/third_party/tiff:tiff",
  ],
)

cc_test(
  name = "frame_test",
  srcs = [ "test/frame_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib" ],
    "//conditions:default": [],
    }) + lcov_linkopts,
  deps = [
    ":frame",
    ":tiff_interface",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ] + select({
    ":win": [ "//bloodflow/third_party/tiff:tiff_dll" ],
    "//conditions:default": [],
  }),
)

cc_library(
  name = "fx3",
  hdrs = [ "fx3.h" ],
  srcs = [ "fx3.cpp" ],
  deps = [
    ":driver",
    ":time",
    "//bloodflow/third_party/cypress-fx3:CyAPI"
  ],
)

cc_test(
  name = "fx3_test",
  srcs = [ "test/fx3_test.cpp" ],
  linkopts = [
    "advapi32.lib",
    "legacy_stdio_definitions.lib",
    "Setupapi.lib",
    "user32.lib",
    "/NODEFAULTLIB:libcmt.lib",
  ],
  deps = [
    ":fx3",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ],
)

cc_library(
  name = "octopus",
  hdrs = [
    "octopus.h",
    "octo_fw.h",
  ],
  srcs = [ "octopus.cpp" ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [ ":driver" ],
)

cc_test(
  name = "octopus_test",
  srcs = [ "test/octopus_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib", "/NODEFAULTLIB:libcmt.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":octopus",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ],
)

cc_library(
  name = "pool",
  srcs = [ "pool.h" ],
)

cc_library(
  name = "rcam",
  hdrs = [
    "rcam.h",
    "rcam_param.h",
    "rcam_fw.h",
  ],
  srcs = [
    "rcam.cpp",
  ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [
    ":circular_buffer",
    ":driver",
    ":execnode",
    ":frame",
    ":time",
  ],
)

cc_test(
  name = "rcam_test",
  srcs = [ "test/rcam_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib", "/NODEFAULTLIB:libcmt.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":rcam",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ] + select({
    ":win": [ "//bloodflow/third_party/tiff:tiff_dll" ],
    "//conditions:default": [],
  }),
)

cc_library(
  name = "stddev",
  hdrs = [ "stddev.h" ],
  srcs = [ "stddev.cpp"],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [
    ":execnode",
    ":frame",
    ":pool",
  ],
)

cc_test(
  name = "stddev_test",
  srcs = [ "test/stddev_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":stddev",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ] + select({
    ":win": [ "//bloodflow/third_party/tiff:tiff_dll" ],
    "//conditions:default": [],
  }),
)

cc_library(
  name = "tiff_interface",
  hdrs = [ "TiffInterface.h" ],
)

cc_library(
  name = "time",
  hdrs = [ "time.h" ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
)

cc_test(
  name = "time_test",
  srcs = [ "test/time_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":time",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ],
)
