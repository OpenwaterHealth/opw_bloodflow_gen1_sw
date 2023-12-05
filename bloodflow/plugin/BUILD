package(default_visibility = ["//visibility:public"])

config_setting(name = "win", constraint_values = ["@platforms//os:windows"])
config_setting(name = "mac", constraint_values = ["@platforms//os:osx"])
config_setting(name = "linux", constraint_values = ["@platforms//os:linux"])

lcov_copts = select({
  ":linux": [ "-fprofile-instr-generate -fcoverage-mapping" ],
  "//conditions:default": [],
})
lcov_linkopts = select({
  ":linux": [ "-fprofile-instr-generate -fcoverage-mapping" ],
  "//conditions:default": [],
})

cc_binary(
  name = "scan.dll",
  linkshared = True,
  srcs = [ "scan.cpp" ],
  linkopts = [
    "advapi32.lib",
    "legacy_stdio_definitions.lib",
    "Setupapi.lib",
    "user32.lib",
    "ws2_32.lib",
    "/NODEFAULTLIB:libcmt.lib",
  ],
  deps = [
    ":CameraManager",
    ":OctopusManager",
    ":VoxelSave",
    "//bloodflow/component:execnode",
    "//bloodflow/component:frame",
    "//bloodflow/component:fx3",
    "//bloodflow/component:octopus",
    "//bloodflow/component:rcam",
    "//bloodflow/component:stddev",
    "//bloodflow/component:time",
    "//bloodflow/third_party/tiff:tiff_dll",
  ],
)

py_test(
  name = "load_tst",
  srcs = [ "load_tst.py" ],
  data = [
    ":scan.dll",
    "//bloodflow/third_party/tiff:tiff_dll",
  ],
)

cc_library(
  name = "CameraManager",
  hdrs = [ "CameraManager.h" ],
  srcs = [ "CameraManager.cpp" ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [
    ":VoxelData",
    ":VoxelSave",
    "//bloodflow/component:rcam",
    "//bloodflow/component:stddev",
    "//bloodflow/third_party/json-develop:json_develop",
  ],
)

cc_test(
  name = 'cameramanager_test',
  srcs = [ "test/CameraManager_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib", "/NODEFAULTLIB:libcmt.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":CameraManager",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ],
)

cc_library(
  name = "OctopusManager",
  hdrs = [ "OctopusManager.h" ],
  srcs = [ "OctopusManager.cpp" ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [
    "//bloodflow/component:octopus",
    "//bloodflow/third_party/json-develop:json_develop",
  ],
)

cc_test(
  name = 'octopusmanager_test',
  srcs = [ "test/OctopusManager_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib", "/NODEFAULTLIB:libcmt.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":OctopusManager",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ],
)

cc_library(
  name = "VoxelData",
  hdrs = [ "VoxelData.h" ],
)

cc_library(
  name = "VoxelSave",
  hdrs = [ "VoxelSave.h" ],
  srcs = [ "VoxelSave.cpp" ],
  copts = lcov_copts,
  linkopts = lcov_linkopts,
  deps = [
    ":VoxelData",
    "//bloodflow/component:execnode",
    "//bloodflow/component:frame",
    "//bloodflow/component:stddev",
  ],
)

cc_test(
  name = "voxelsave_test",
  srcs = [ "test/voxelsave_test.cpp" ],
  copts = lcov_copts,
  linkopts = select({
    ":win": [ "advapi32.lib", "user32.lib" ],
    "//conditions:default": [],
  }) + lcov_linkopts,
  deps = [
    ":VoxelSave",
    "//googletest:gtest",
    "//googletest:gtest_main",
  ],
)
