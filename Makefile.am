warning_flags = -Wall -Wextra -Wshadow -Wpedantic
common_cflags = -g3 -O0 $(MFLAGS) $(STACKREALIGN) $(warning_flags)
AM_CXXFLAGS = -g3 -std=c++11 $(common_cflags)

AM_CPPFLAGS = $(VapourSynth_CFLAGS)

lib_LTLIBRARIES = libhqdn3d.la

libhqdn3d_la_SOURCES = src/hqdn3d.cpp

libhqdn3d_la_LDFLAGS = -no-undefined -avoid-version $(PLUGINLDFLAGS)
