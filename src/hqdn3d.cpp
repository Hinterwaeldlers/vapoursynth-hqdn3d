/*
    HQDN3D 1.00 for Vapoursynth

    Copyright (C) 2003 Daniel Moreno <comac@comac.darktech.org>
    Avisynth port (C) 2005 Loren Merritt <lorenm@u.washington.edu>
    Vapoursynth port (C) 2017 Martin Güthle  <mguethle@xunit.de>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <algorithm>

#include <VapourSynth.h>
#include <VSHelper.h>

typedef struct Hqdn3dData {
    VSNodeRef *clip;
    const VSVideoInfo *vi;

    double lumSpac;
    double chromSpac;
    double lumTmp;
    double chromTmp;
    int restartLap;

    int coefs[4][512*16];
    unsigned int *prevFrame[3] = { nullptr, nullptr, nullptr };
    unsigned int *prevLine[3] = { nullptr, nullptr, nullptr };
    bool process[3];
    int last_frame = -1; // the last frame returned by hqdn3d
} Hqdn3dData;

static void
VS_CC hqdn3dFree(void *instanceData, VSCore *core, const VSAPI *vsapi) {
    (void)core;

    Hqdn3dData *d = (Hqdn3dData *)instanceData;

    for (int p = 0; p < d->vi->format->numPlanes; p++) {
        free(d->prevLine[p]);
        free(d->prevFrame[p]);
    }

    vsapi->freeNode(d->clip);
    free(d);
}

static inline unsigned int
LowPassMul(unsigned int pMul, unsigned int cMul, const int* coef){
    static const unsigned int ROUND_CONVOLUTION = 0x10007FF;
    static const unsigned int SHIFT_CONVOLUTION = 12;
    int d = (static_cast<int>(pMul - cMul) + ROUND_CONVOLUTION) >> SHIFT_CONVOLUTION;
    return cMul + coef[d];
}

static void
deNoise(
      const uint8_t * srcPlane
    , unsigned int * prevPlane
    , unsigned int * prevLine
    , uint8_t * tarPlane
    , const int frameWidth
    , const int frameHeight
    , const int srcStride
    , const int tarStride
    , const int *coefsHorizontal
    , const int *coefsVertical
    , const int *coefsTemporal
    , const bool isFirstFrame
) {
    static const unsigned int ROUND_LINE  = 0x1000007F;
    static const unsigned int SHIFT_LINE  = 8;
    static const unsigned int ROUND_PIXEL = 0x10007FFF;
    static const unsigned int SHIFT_PIXEL = 16;

    for (int row = 0; row < frameHeight; ++row) {
        /* gcc assume prevPixel might be used in an uninitialized way, but
         * it's not. So feel free to use any other value
         */
        unsigned int prevPixel = 0;
        for (int col = 0; col < frameWidth; ++col) {
            // Correlate current pixel with previous pixel
            prevPixel = col == 0
                ? srcPlane[row * srcStride + col] << SHIFT_PIXEL
                : LowPassMul(
                      prevPixel
                    , srcPlane[row * srcStride + col] << SHIFT_PIXEL
                    , coefsHorizontal
                );
            // Correlate previous line with previous pixel
            prevLine[col] = row == 0
                ? prevPixel
                : LowPassMul(
                      prevLine[col]
                    , prevPixel
                    , coefsVertical
                );
            unsigned int resPix;
            if (isFirstFrame) {
                resPix = prevLine[col];
            } else {
                // Correlate vertical result with previous result frame pixel
                resPix = LowPassMul(
                      prevPlane[row * frameWidth + col] << SHIFT_LINE
                    , prevLine[col]
                    , coefsTemporal
                );
            }

            prevPlane[row * frameWidth + col]
                = ((resPix + ROUND_LINE) >> SHIFT_LINE) & 0xFFFF;
            if (tarPlane)
                tarPlane[row * tarStride + col]
                    = (resPix + ROUND_PIXEL) >> SHIFT_PIXEL;
        }
    }
}

static void filterFrame(
      const VSFrameRef *srcFrame
    , VSFrameRef *newFrame
    , const bool isFirstFrame
    , Hqdn3dData *usrData
    , const VSAPI *vsapi
) {
    const VSFormat *srcFrameFmt = vsapi->getFrameFormat(srcFrame);

    for (int plane = 0; plane < srcFrameFmt->numPlanes; plane++) {
        if (!usrData->process[plane])
            continue;

        deNoise(
              vsapi->getReadPtr(srcFrame, plane)
            , usrData->prevFrame[plane]
            , usrData->prevLine[plane]
            , newFrame ? vsapi->getWritePtr(newFrame, plane) : nullptr
            , vsapi->getFrameWidth(srcFrame, plane)
            , vsapi->getFrameHeight(srcFrame, plane)
            , vsapi->getStride(srcFrame, plane)
            , newFrame ? vsapi->getStride(newFrame, plane) : 0
            , usrData->coefs[plane == 0 ? 0 : 2] // Y or U/V
            , usrData->coefs[plane == 0 ? 0 : 2] // Y or U/V
            , usrData->coefs[plane == 0 ? 1 : 3] // Y or U/V
            , isFirstFrame
        );
    }

}

static const VSFrameRef *VS_CC hqdn3dGetFrame(
      int n
    , int activationReason
    , void **instanceData
    , void **frameData
    , VSFrameContext *frameCtx
    , VSCore *core
    , const VSAPI *vsapi
) {
    (void)frameData;

    // Get the user data
    Hqdn3dData * usrData = reinterpret_cast<Hqdn3dData *>(*instanceData);

    if (activationReason == arInitial) {
        // if we skip some frames, filter the gap anyway
        if (n > usrData->last_frame + 1 &&
            n - usrData->last_frame <= usrData->restartLap + 1 &&
            usrData->last_frame >= 0) {

            for (int i = usrData->last_frame + 1; i < n; i++) {
                vsapi->requestFrameFilter(i, usrData->clip, frameCtx);
            }
        // if processing out of sequence, filter several previous frames to minimize seeking problems
        } else if (n != usrData->last_frame + 1) {
            int sn = std::max(0, n - usrData->restartLap);

            for (int i = sn + 1; i < n; i++)
                vsapi->requestFrameFilter(i, usrData->clip, frameCtx);
        }

        vsapi->requestFrameFilter(n, usrData->clip, frameCtx);

        return nullptr;
    }
    if (activationReason != arAllFramesReady) {
        return nullptr;
    }

    // if we skip some frames, filter the gap anyway
    if (n > usrData->last_frame + 1 &&
        n - usrData->last_frame <= usrData->restartLap + 1 &&
        usrData->last_frame >= 0) {

        for (int i = usrData->last_frame + 1; i < n; i++) {
            const VSFrameRef *f = vsapi->getFrameFilter(i, usrData->clip, frameCtx);

            filterFrame(f, nullptr, false, usrData, vsapi);

            vsapi->freeFrame(f);
        }
    // if processing out of sequence, filter several previous frames to minimize seeking problems
    } else if (n != usrData->last_frame + 1) {
        int sn = std::max(0, n - usrData->restartLap);

        for (int i = sn + 1; i < n; i++) {
            const VSFrameRef *f = vsapi->getFrameFilter(i, usrData->clip, frameCtx);

            filterFrame(f, nullptr, i == sn + 1, usrData, vsapi);

            vsapi->freeFrame(f);
        }
    }


    // Get current frame
    const VSFrameRef * srcFrame =
        vsapi->getFrameFilter(n, usrData->clip, frameCtx);


    // Create target frame
    const VSFrameRef *plane_src[3] = {
        usrData->process[0] ? nullptr : srcFrame,
        usrData->process[1] ? nullptr : srcFrame,
        usrData->process[2] ? nullptr : srcFrame
    };
    int planes[3] = { 0, 1, 2 };

    VSFrameRef *newFrame = vsapi->newVideoFrame2(
          usrData->vi->format
        , usrData->vi->width
        , usrData->vi->height
        , plane_src
        , planes
        , srcFrame
        , core
    );

    filterFrame(srcFrame, newFrame, n == 0, usrData, vsapi);

    vsapi->freeFrame(srcFrame);

    usrData->last_frame = n;

    return newFrame;
}


static void VS_CC hqdn3dInit(
      VSMap *in
    , VSMap *out
    , void **instanceData
    , VSNode *node
    , VSCore *core
    , const VSAPI *vsapi
) {
    (void)in;
    (void)out;
    (void)core;

    Hqdn3dData *d = (Hqdn3dData *) * instanceData;

    vsapi->setVideoInfo(d->vi, 1, node);
}

// Create the plugin
static void VS_CC hqdn3dCreate(
      const VSMap *in
    , VSMap *out
    , void *userData
    , VSCore *core
    , const VSAPI *vsapi
) {
    (void)userData;

    Hqdn3dData d;
    int err;

    d.lumSpac    = vsapi->propGetFloat(in, "lum_spac",    0, &err);
    if (err) {
        d.lumSpac = 4.0;
    } else if (d.lumSpac < 0 || d.lumSpac > 255) {
        vsapi->setError(out, "Hqdn3d: lum_spac must be between 0 and 255 (inclusive).");
        return;
    }

    d.chromSpac  = vsapi->propGetFloat(in, "chrom_spac",  0, &err);
    if (err) {
        d.chromSpac = .75 * d.lumSpac;
    } else if (d.chromSpac < 0 || d.chromSpac > 255) {
        vsapi->setError(out, "Hqdn3d: chrom_spac must be between 0 and 255 (inclusive).");
        return;
    }

    d.lumTmp     = vsapi->propGetFloat(in, "lum_tmp",     0, &err);
    if (err) {
        d.lumTmp = 1.5 * d.lumSpac;
    } else if (d.lumTmp < 0 || d.lumTmp > 255) {
        vsapi->setError(out, "Hqdn3d: lum_tmp must be between 0 and 255 (inclusive).");
        return;
    }

    d.chromTmp   = vsapi->propGetFloat(in, "chrom_tmp",   0, &err);
    if (err) {
        d.chromTmp = (d.lumSpac == 0) ? d.chromSpac * 1.5
                                      : d.lumTmp * d.chromSpac / d.lumSpac;
    } else if (d.chromTmp < 0 || d.chromTmp > 255) {
        vsapi->setError(out, "Hqdn3d: chrom_tmp must be between 0 and 255 (inclusive).");
        return;
    }

    d.restartLap = int64ToIntS(vsapi->propGetInt(in, "restart_lap", 0, &err));
    if (err)
        d.restartLap = std::max(2
            , static_cast<int>(1 + std::max(d.lumTmp, d.chromTmp)));


    d.clip = vsapi->propGetNode(in, "clip", 0, NULL);
    d.vi = vsapi->getVideoInfo(d.clip);

    if (!d.vi->format ||
        d.vi->format->colorFamily == cmRGB ||
        d.vi->format->bitsPerSample > 8 ||
        d.vi->width == 0 ||
        d.vi->height == 0) {

        vsapi->setError(out, "Hqdn3d: input clip must be 8 bit, not RGB, and it must have constant format and dimensions.");
        vsapi->freeNode(d.clip);
        return;
    }


    d.lumSpac   = std::min(254.9, d.lumSpac);
    d.chromSpac = std::min(254.9, d.chromSpac);
    d.lumTmp    = std::min(254.9, d.lumTmp);
    d.chromTmp  = std::min(254.9, d.chromTmp);

    // Calculate the coefficients
    for (auto const cc : {
          std::make_pair(0, d.lumSpac)
        , std::make_pair(1, d.lumTmp)
        , std::make_pair(2, d.chromSpac)
        , std::make_pair(3, d.chromTmp)
    } ) {
        const double gamma = std::log(0.25) / std::log(1.0 - cc.second / 255.0 - 0.00001);
        for (int i = -255 * 16; i < 256 * 16; ++i) {
            const double simil = 1.0 - std::abs(i) / (16*255.0);
            const double c = std::pow(simil, gamma) * 65536.0 * i / 16.0;
            d.coefs[cc.first][16*256+i]
                = static_cast<int>(c < 0 ? c - 0.5 : c + 0.5);
        }
    }

    // According to the documentation, 0 strength means no processing.
    d.process[0] = d.lumSpac != 0 || d.lumTmp != 0;
    d.process[1] = d.process[2] = d.chromSpac != 0 || d.chromTmp != 0;


    for (int p = 0; p < d.vi->format->numPlanes; p++) {
        int width = d.vi->width;
        int height = d.vi->height;

        if (p) {
            width >>= d.vi->format->subSamplingW;
            height >>= d.vi->format->subSamplingH;
        }

        d.prevFrame[p] = (unsigned int *)malloc(width * height * sizeof(unsigned int));
        d.prevLine[p] = (unsigned int *)malloc(width * sizeof(unsigned int));
    }


    Hqdn3dData * data = (Hqdn3dData *)malloc(sizeof(Hqdn3dData));
    *data = d;

    vsapi->createFilter(
          in
        , out
        , "Hqdn3d"
        , hqdn3dInit
        , hqdn3dGetFrame
        , hqdn3dFree
        , fmSerial
        , nfMakeLinear
        , data
        , core
    );
}

VS_EXTERNAL_API(void) VapourSynthPluginInit(
      VSConfigPlugin configFunc
    , VSRegisterFunction registerFunc
    , VSPlugin *plugin
) {
    configFunc(
          "com.vapoursynth.hqdn3d"
        , "hqdn3d"
        , "HQDn3D port as used in avisynth/mplayer"
        , VAPOURSYNTH_API_VERSION
        , 1
        , plugin
    );
    registerFunc(
          "Hqdn3d"
        , "clip:clip;"
          "lum_spac:float:opt;"
          "chrom_spac:float:opt;"
          "lum_tmp:float:opt;"
          "chrom_tmp:float:opt;"
          "restart_lap:int:opt;"
        , hqdn3dCreate
        , 0
        , plugin
    );
}
