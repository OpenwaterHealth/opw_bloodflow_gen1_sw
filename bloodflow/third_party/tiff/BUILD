package(default_visibility = ["//visibility:public"])

config_setting(name = "win", constraint_values = ["@platforms//os:windows"])
config_setting(name = "linux", constraint_values = ["@platforms//os:linux"])

cc_library(
  name = "tiff",
  hdrs = [
    "inc/tiff.h",
    "inc/tiffconf.h",
    "inc/tiffio.h",
    "inc/tiffvers.h",
  ],
  srcs = select({
    ":win": [ "win/lib/tiff.lib" ],
    ":linux": [ "linux/lib/libtiff.so.5", "linux/lib/libjbig.so.0" ],
    "//conditions:default": []
  }),
)

cc_import(
  name = "tiff_dll",
  shared_library = "win/dll/tiff.dll",
)
