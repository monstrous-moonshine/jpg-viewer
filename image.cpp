#include <iostream>
#include "image.h"
#include <cmath>

Uint8* copy(Uint8* dest, Uint8* source, int width, int height) {
  std::cerr << "Copying...";
  for (int i = 0; i < width * height * 4; i++) {
    dest[i] = source[i];
  }
  std::cerr << "done.\n";
  
  return dest;
}

float linear_interp(float x, float p0, float p1) {
    return (1 - x) * p0 + x * p1;
}

float cubic_interp(float x, float p0, float p1, float p2, float p3) {
    return p1 + 0.5 * x * (p2 - p0 + x * (2 * p0 - 5 * p1 + 4 * p2 - p3 + x * (
                                             -p0 + 3 * (p1 - p2) + p3)));
}

// bilinear or bicubic interpolation
Uint8* rescale(Uint8* source, int sourcew, int sourceh, int destw, int desth) {
  Uint8* dest = new Uint8[destw * desth * 4];
  if (sourcew == destw) {
    return copy(dest, source, sourcew, sourceh);
  }
  
  // this is dangerous, and only works because the two
  // ratios are slightly different, and the lesser of
  // the two makes the last index miss the mark for a little
  // to be safe, subtract 4, and then you can afford to
  // use max instead of min
  float scale = fmin((sourcew - 3) / float(destw - 1),
                    (sourceh - 3) / float(desth - 1)
                    );
  // resample
  std::cerr << "Interpolating...";
  for (int i = 0; i < desth; i++) {
    for (int j = 0; j < destw; j++) {
      for (int k = 0; k < 3; k++) {
        int m = 1 + i * scale;
        int n = 1 + j * scale;
        // interpolate
        float y = i * scale - (m - 1);
        float x = j * scale - (n - 1);
#define SRC(m, n) (source[((m) * sourcew + n) * 4 + k])
        /*
#define SRC2(m, n) SRC(m, n), SRC(m, n + 1)
        // bilinear
        float r = linear_interp(x, SRC2(m    , n));
        float s = linear_interp(x, SRC2(m + 1, n));
        float v = linear_interp(y, r, s);
#undef SRC2
        */
#define SRC4(m, n) SRC(m, n - 1), SRC(m, n), SRC(m, n + 1), SRC(m, n + 2)
        // bicubic
        float r = cubic_interp(x, SRC4(m - 1, n));
        float s = cubic_interp(x, SRC4(m    , n));
        float t = cubic_interp(x, SRC4(m + 1, n));
        float u = cubic_interp(x, SRC4(m + 2, n));
        float v = cubic_interp(y, r, s, t, u);
        v = v > 255 ? 255 : (v < 0 ? 0 : v);
#undef SRC4
#undef SRC
#define DST(i, j, k) dest[((i) * destw + j) * 4 + k]
        DST(i, j, k) = v;
      }
      
      // copy the alpha channel
      DST(i, j, 3) = 255;
#undef DST
    }
  }
  std::cerr << "done.\n";
  
  return dest;
}
