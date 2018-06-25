Description
===========

High Quality DeNoise 3D performs a 3-way low-pass filter, which can
completely remove high-frequency noise while minimizing blending
artifacts.

Due to the nature of the filter, there is no possibility of
multithreading. Scripts where this filter is used most likely won't
keep a multi-core CPU fully busy. If this is a problem, the video
should be divided into a few chunks at scene changes and filtered with
several instances of VapourSynth.

This is a port of the Hqdn3d Avisynth plugin version 0.11, which is in
turn based on the filter with the same name from mplayer.


Usage
=====
::

    hqdn3d.Hqdn3d(clip clip, [float lum_spac, float chrom_spac, float lum_tmp, float chrom_tmp, int restart_lap])


Parameters:
    *clip*
        A clip to process. It must have constant format and dimensions
        and it must be 8 bit and not RGB.

    *lum_spac*
        Luma spatial filter strength. Must be between 0 and 255.

        Default: 4.0.

    *chrom_spac*
        Chroma spatial filter strength. Must be between 0 and 255.

        Default: lum_spac * 0.75.

    *lum_tmp*
        Luma temporal filter strength. Must be between 0 and 255.

        Default: lum_spac * 1.5.

    *chrom_tmp*
        Chroma temporal filter strength. Must be between 0 and 255.

        Default when lum_spac is 0: chrom_spac * 1.5.

        Default when lum_spac is not 0: lum_tmp * chrom_spac / lum_spac.

    *restart_lap*
        Whenever a frame is requested out of order, restart filtering
        this many frames before. While seeking still slightly affects
        the content of the frames returned, this should reduce the
        disturbance to an unnoticable level. 

        Default: max(2, int(1 + max(lum_tmp, chrom_tmp)))


Compilation
===========

::

    ./autogen.sh
    ./configure
    make


License
=======

GNU GPL v2, like the Avisynth plugin.
