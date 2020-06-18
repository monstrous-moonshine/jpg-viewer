#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <memory>
#include <limits>
#include <cmath>
#include "jpeg.h"

std::pair<int, unsigned char*> readFile(const char* path) {
  using std::ios;
  std::ifstream file(path, ios::in | ios::binary | ios::ate);
  
  if (file.is_open()) {
    int fileSize = file.tellg();
    unsigned char* buffer = new unsigned char[fileSize];
    file.seekg(0, ios::beg);
    file.read((char*) buffer, fileSize);
    return std::make_pair(fileSize, buffer);
  } else {
    std::cerr << "Error opening file.\n";
    exit(EXIT_FAILURE);
  }
}

void image_t::print_markers() {
  std::cerr << "Position Marker\n";
  for (int i = 0; i < fileSize; i++) {
    if (buffer[i] == 0xff && buffer[i+1] != 0 && buffer[i+1] != 0xff) {
      std::cerr << std::setw(8) << i << " ";
      std::cerr << std::setw(2);
      std::cerr << std::hex << std::setfill('0');
      std::cerr << int(buffer[i+1]);
      std::cerr << std::dec << std::setfill(' ') << "\n";
    }
  }
  std::cerr << std::setw(0) << "\n"; // restore default width
}

image_t::image_t(
  const char* filePath,
  int width  = std::numeric_limits<int>::max(),
  int height = std::numeric_limits<int>::max()
                ) {
  std::pair<int, unsigned char*> file = readFile(filePath);
  fileSize = file.first;
  buffer = file.second;
  print_markers();
  dispWidth = width;
  dispHeight = height;
  read_jpeg();
}

image_t::~image_t() {
  for (int i = 0; i < 4; i++) {
    delete htables[i];
  }
  for (int i = 0; i < 2; i++) {
    delete qtables[i];
  }
  delete pixels;
}

void image_t::read_jpeg() {
  read_marker(SOI);
  read_tables();
  read_marker(SOF);
  read_marker(EOI);
  dequantize();
  unzigzag();
  idct();
  color_convert();
  pixel_rearrange();
}

void image_t::read_marker(int marker_type) {
  int i = this->pos;
  
  // skip over fill bytes, and enforce correct header
  while (buffer[i] == 0xff && buffer[i+1] == 0xff) i++;
  if (buffer[i] != 0xff || buffer[i+1] != marker_type) {
    std::cerr << "Marker ";
    std::cerr << std::hex << marker_type;
    std::cerr << " not found at position ";
    std::cerr << std::dec << i << ".\n";
    exit(EXIT_FAILURE);
  }
  i += 2, this->pos = i;
  
  // check for markers without segments
  if (marker_type == SOI) return;
  if (marker_type == EOI) {
    delete buffer;
    return;
  }
  
  // read marker segment
  int len = buffer[i] * 256 + buffer[i+1];
  i += 2, this->pos = i;
  switch (marker_type) {
  case SOF:
    read_frame();
    break;
  case SOS:
    read_scan();
    break;
  case DHT:
    read_hufftable(len - 2);
    break;
  case DQT:
    read_quantable(len - 2);
    break;
  case DRI:
    read_rst();
    break;
  case COM:
    this->pos += len - 2;
    break;
  default:
    if (marker_type >= APP && marker_type < JPG) {
      this->pos += len - 2;
    } else {
      std::cerr << "Unrecognized marker at position " << i << ".\n";
      exit(EXIT_FAILURE);
    }
  }
}

void image_t::read_tables() {
  int i = this->pos;
  
  // skip over fill bytes and check for marker
  while (buffer[i] == 0xff && buffer[i+1] == 0xff) i++;
  if (buffer[i] != 0xff) {
    std::cerr << "Marker not found at position " << i << ".\n";
    exit(EXIT_FAILURE);
  }
  this->pos = i;
  
  // process tables, ignore others
  unsigned char marker_type = buffer[i+1];
  switch (marker_type) {
  case DHT:
  case DQT:
  case DRI:
  case COM:
    read_marker(marker_type);
    read_tables();
    break;
  case DAC:
    std::cerr << "Arithmetic coding not supported.\n";
    exit(EXIT_FAILURE);
  default:
    if (marker_type >= APP && marker_type < JPG) {
      // application data, ignore
      read_marker(marker_type);
      read_tables();
    }
  }
}

void image_t::read_hufftable(int len) {
  int i = this->pos;
  
  while (i - this->pos < len) {
    int Tc = buffer[i] / 16;
    int Th = buffer[i] % 16;
    i++;
    int bits[17];
    int total = 0;
    for (int nbits = 1; nbits <= 16; nbits++) {
      bits[nbits] = buffer[i++];
      total += bits[nbits];
    }
    
    std::unique_ptr<int[]> huffval(new int[total]);
    for (int j = 0; j < total; j++) {
      huffval[j] = buffer[i++];
    }
    
    this->htables[2*Tc+Th] = make_hufftable(bits, huffval.get());
  }
  
  this->pos = i;
}

void image_t::read_quantable(int len) {
  int i = this->pos;
  assert(len % 65 == 0); // 1 byte for Pq and Tq, 64 for data
  
  for (int table = 0; table < len / 65; table++) {
    int Pq = buffer[i] / 16;
    int Tq = buffer[i] % 16;
    assert(Pq == 0); // 8 bit precision
    assert(Tq < 2);  // 2 tables, luma and chroma
    i++;
    this->qtables[Tq] = new int[64];
    for (int j = 0; j < 64; j++) {
      this->qtables[Tq][j] = buffer[i++];
    }
  }
  
  this->pos = i;
}

void image_t::read_rst() {
  int i = this->pos;
  std::cerr << "Enabling restart interval...";
  this->rst = buffer[i] * 256 + buffer[i + 1];
  std::cerr << "restart interval is " << rst << ".\n";
  this->pos += 2;
}

// adjust this to choose between accuracy and speed
int calculateQuality(int srcw, int srch, int dstw, int dsth) {
  float scale = fmax(srcw  / float(dstw), srch / float(dsth));
  if (scale < 2) {
    return 8;
  } else if (scale < 4) {
    return 4;
  } else if (scale < 8) {
    return 2;
  } else {
    return 1;
  }
}

void image_t::read_frame() {
  int i = this->pos;
  int P = buffer[i++];
  assert(P == 8); // precision
  this->height = buffer[i] * 256 + buffer[i + 1];
  i += 2;
  this->width = buffer[i] * 256 + buffer[i + 1];
  i += 2;
  assert(height > 0); // do not support DNL format
  int Nf = buffer[i++];
  assert(Nf == 3);    // JFIF Y, Cb, Cr components
  for (int j = 0; j < Nf; j++, i += 3) {
    Ci[j] = buffer[i+0];
    Hi[j] = buffer[i+1] / 16;
    Vi[j] = buffer[i+1] % 16;
    Tq[j] = buffer[i+2];
    assert(Tq[j] < 2);
    hmax = Hi[j] > hmax ? Hi[j] : hmax;
    vmax = Vi[j] > vmax ? Vi[j] : vmax;
    size_mcu += Hi[j] * Vi[j];
    std::cerr << "Component " << Ci[j] << ": ";
    std::cerr << Hi[j] << "x" << Vi[j] << "\n";
  }
  quality = calculateQuality(width, height, dispWidth, dispHeight);
  assert(size_mcu <= 10); // JPEG specification
  this->nmcus = ceil(width/(8.0*hmax)) * ceil(height/(8.0*vmax));
  this->rst = rst > 0 ? rst : nmcus;
  this->data   = new int[nmcus * size_mcu * 64];
  this->padded = new unsigned char[nmcus * hmax * vmax * 4 * 64];
  this->pos = i;
  read_tables();
  read_marker(SOS);
}

char nextbit(unsigned char* buffer, int& index, bool restart = false) {
  static unsigned char byte;
  static int count = -1;
  
  // if at the end of entropy coded segment
  // reset count to empty
  if (restart) {
    count = -1;
    return '\0';
  }
  
  // if run out of bits, fetch next byte
  if (count < 0) {
    byte = buffer[index++];
    
    // check for markers
    if (byte == 0xff) {
      if (buffer[index] != 0) {
        std::cerr << "Unexpected marker in entropy coded segment: ";
        std::cerr << std::hex << int(buffer[index]);
        std::cerr << " at position ";
        std::cerr << std::dec << index << ".\n";
        exit(EXIT_FAILURE);
      } else {
        // this is a stuffing zero byte, ignore
        index++;
      }
    }
    
    // reset count to full
    count = 7;
  }
  
  // update count and return bit
  return (byte >> count--) & 1;
}

int readEntry(int nbits, unsigned char* buffer, int& index) {
  if (nbits == 0) {
    return 0;
  }
  int entry = 0;
  for (int i = 0; i < nbits; i++) {
    entry <<= 1;
    entry += nextbit(buffer, index);
  }
  if (entry < (1 << (nbits - 1))) {
    entry = -((1 << nbits) - 1 - entry);
  }
  return entry;
}

int readValue(unsigned char* buffer, huffnode* htable, int& i) {
  int idx = 0;
  int value = htable[idx].value;
  while (value == NONTERMINAL) {
    char c = nextbit(buffer, i);
    idx = htable[idx].child[c];
    value = htable[idx].value;
  }
  return value;
}

void image_t::read_scan() {
  int i = this->pos;
  int Ns = buffer[i++];
  assert(Ns == 3);      // Ns = 3 for JFIF
  int Cs[3], Td[3], Ta[3];
  for (int j = 0; j < 3; j++, i++) {
    Cs[j] = buffer[i++];
    Td[j] = buffer[i] / 16;
    Ta[j] = buffer[i] % 16;
  }
  i += 3; // skip 00 3f 00
  
  // array of previously decoded DC values
  int prev[3];
  
  std::cerr << "Decoding...";
  for (int rstint = 0, mcu = 0; mcu < nmcus; rstint++) { // rstint = restart interval
    
    // init prev array
    for (int j = 0; j < 3; j++) {
      prev[j] = 0;
    }
    
    // read MCUs
    for (int rstmcu = 0 ; rstmcu < this->rst; rstmcu++, mcu++) {
      // read components
      for (int comp = 0, cumdu = 0; comp < 3; comp++) { // cumdu = cumulative du
        // read data units
        for (int du = 0; du < Hi[comp] * Vi[comp]; du++, cumdu++) {
#define DU (mcu * size_mcu + cumdu)
          int nElems = 0;
          // read DC value
          {
            int value = readValue(buffer, htables[Td[comp]], i);
            data[DU * 64] = readEntry(value, buffer, i) + prev[comp];
            prev[comp] = data[DU * 64 + nElems++];
          }
          // read AC values
          while (nElems < 64) {
            int value = readValue(buffer, htables[2+Ta[comp]], i);
            if (value == 0) {
              // EOB
              while (nElems < 64) {
                data[DU * 64 + nElems++] = 0;
              }
            } else if (value == 0xf0) {
              // ZRL
              for (int j = 0; j < 16; j++) {
                data[DU * 64 + nElems++] = 0;
              }
            } else {
              int nzero = value / 16;
              int nbits = value % 16;
              for (int j = 0; j < nzero; j++) {
                data[DU * 64 + nElems++] = 0;
              }
              data[DU * 64 + nElems++] = readEntry(nbits, buffer, i);
            }
          }
#undef DU
        }
      }
    }
    
    // read RST marker if not last entropy coded segment
    if (mcu < nmcus) {
      if (buffer[i] != 0xff || buffer[i + 1] != RST + rstint % 8) {
        std::cerr << "Error reading restart marker at position ";
        std::cerr << i << ".\n";
        exit(EXIT_FAILURE);
      } else {
        nextbit(buffer, i, true); // restart
        i += 2;
      }
    }
  }
  std::cerr << "done.\n";
  
  this->pos = i;
}

void image_t::dequantize() {
  for (int mcu = 0; mcu < nmcus; mcu++) {
    for (int comp = 0, cumdu = 0; comp < 3; comp++) {
      for (int du = 0; du < Hi[comp] * Vi[comp]; du++, cumdu++) {
#define DU (mcu*size_mcu+cumdu)
        for (int j = 0; j < 64; j++) {
          data[DU * 64 + j] *= qtables[Tq[comp]][j];
        }
#undef DU
      }
    }
  }
}

void image_t::unzigzag() {
  int temparray[64];
  
  int zigzagx[64] = {
    0,
    1, 0,
    0, 1, 2,
    3, 2, 1, 0,
    0, 1, 2, 3, 4,
    5, 4, 3, 2, 1, 0,
    0, 1, 2, 3, 4, 5, 6,
    7, 6, 5, 4, 3, 2, 1, 0,
    1, 2, 3, 4, 5, 6, 7,
    7, 6, 5, 4, 3, 2,
    3, 4, 5, 6, 7,
    7, 6, 5, 4,
    5, 6, 7,
    7, 6,
    7
  };
  
  int zigzagy[64] = {
    0,
    0, 1,
    2, 1, 0,
    0, 1, 2, 3,
    4, 3, 2, 1, 0,
    0, 1, 2, 3, 4, 5,
    6, 5, 4, 3, 2, 1, 0,
    0, 1, 2, 3, 4, 5, 6, 7,
    7, 6, 5, 4, 3, 2, 1,
    2, 3, 4, 5, 6, 7,
    7, 6, 5, 4, 3,
    4, 5, 6, 7,
    7, 6, 5,
    6, 7,
    7
  };
  
  // loop over data units
  for (int du = 0; du < nmcus * size_mcu; du++) {
    // copy data into temporary array
    for (int j = 0; j < 64; j++) {
      temparray[j] = data[du*64+j];
    }
    
    // copy back in correct order
    for (int j = 0; j < 64; j++) {
      int x = zigzagx[j];
      int y = zigzagy[j];
      data[du*64+8*y+x] = temparray[j];
    }
  }
}

#define round(x) (int) (x + 0.5)

void image_t::idct() {
  float result[8][8];
  float transmat[8][8];
  float tempmat[8][8];
  
  // compute transformation (Fourier) matrix
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 8; j++) {
      transmat[i][j] = cos(i*(2*j+1)*M_PI/16) * (i > 0 ? 1 : 1 / sqrt(2));
    }
  }
  
  // transform each data unit
  std::cerr << "Inverse discrete cosine transform...";
  for (int mcu = 0, cumdu = 0; mcu < nmcus; mcu++) {
    for (int comp = 0; comp < 3; comp++) {
      for (int du = 0; du < Hi[comp] * Vi[comp]; du++, cumdu++) {
        int hq = fmin(quality * hmax / Hi[comp], 8);
        int vq = fmin(quality * vmax / Vi[comp], 8);
        
        // compute idct
        for (int p = 0; p < 8; p += 8 / vq) {
          for (int q = 0; q < quality; q++) {
            tempmat[p][q] = 0;
            for (int r = 0; r < 8; r++) {
              tempmat[p][q] += transmat[r][p] * data[cumdu*64+8*r+q];
            }
          }
        }
        
        for (int p = 0; p < 8; p += 8 / vq) {
          for (int q = 0; q < 8; q += 8 / hq) {
            result[p][q] = 0;
            for (int r = 0; r < quality; r++) {
              result[p][q] += tempmat[p][r] * transmat[r][q];
            }
          }
        }
        
        // scale, shift, clip
        for (int p = 0; p < 8; p += 8 / vq) {
          for (int q = 0; q < 8; q += 8 / hq) {
            result[p][q] /= 4;
            result[p][q] += 128;
            data[cumdu*64+8*p+q] = round(result[p][q]);
          }
        }
      }
    }
  }
  std::cerr << "done.\n";
}

#define clip(x) x = x < 0 ? 0 : (x > 255 ? 255 : x)

void image_t::color_convert() {
  // index of each component's data units
  int comp[3];
  comp[0] = 0;
  comp[1] = Hi[0] * Vi[0] + comp[0];
  comp[2] = Hi[1] * Vi[1] + comp[1];
  
  // horizontal and vertical subsampling factors
  int hsf[3], vsf[3];
  for (int i = 0; i < 3; i++) {
    hsf[i] = hmax / Hi[i];
    vsf[i] = vmax / Vi[i];
  }
  
#define SRC(i) (mcu * size_mcu + comp[i]  \
                    + vert/vsf[i] * Hi[i] \
                    + horz/hsf[i])
#define IND(i) (8 * ((8/vsf[i]*vert) % 8 + l/vsf[i]) \
                  +  (8/hsf[i]*horz) % 8 + m/hsf[i])
#define DST    ((mcu * vmax + vert) * hmax + horz)

  // convert
  std::cerr << "Color conversion...";
  for (int mcu = 0; mcu < nmcus; mcu++) {
    for (int vert = 0; vert < vmax; vert++) {
      for (int horz = 0; horz < hmax; horz++) {
        for (int i = 0; i < 8; i += 8 / quality) {
          for (int j = 0; j < 8; j += 8 / quality) {
            int k = 8 * i + j;
            int l = k / 8;
            int m = k % 8;
            
            // read data
            int  y = data[SRC(0)*64+IND(0)];
            int cb = data[SRC(1)*64+IND(1)] - 128;
            int cr = data[SRC(2)*64+IND(2)] - 128;
            
            // convert
            float r = y                 +   1.402f * cr;
            float g = y - 0.34414f * cb - 0.71414f * cr;
            float b = y +   1.772f * cb;
            
            // clip values to within range
            clip(r); clip(g); clip(b);
            
            // transfer to padded array
            padded[(DST*64+k)*4+0] = round(r);
            padded[(DST*64+k)*4+1] = round(g);
            padded[(DST*64+k)*4+2] = round(b);
            padded[(DST*64+k)*4+3] = 255;
          }
        }
      }
    }
  }
  std::cerr << "done.\n";
  
  delete this->data;
#undef SRC
#undef IND
#undef DST
}

void image_t::pixel_rearrange() {
#define BLK ((ymcu * nxmcus + xmcu) * \
  hmax * vmax + vblk * hmax + hblk)
#define q quality
  // adjusted height and width
  aWidth = ceil(width * q / 8.0);
  aHeight = ceil(height * q / 8.0);
  this->pixels = new unsigned char[aWidth * aHeight * 4];
  int nxmcus = ceil(width /(8.0*hmax));
  int nymcus = ceil(height/(8.0*vmax));
  
  // rearrange
  std::cerr << "Pixel rearrangement...";
  for (int ymcu = 0; ymcu < nymcus; ymcu++) {
    for (int xmcu = 0; xmcu < nxmcus; xmcu++) {
      for (int vblk = 0; vblk < vmax; vblk++) {
        for (int hblk = 0; hblk < hmax; hblk++) {
          for (int y = 0; y < 8; y += 8 / q) {
            for (int x = 0; x < 8; x += 8 / q) {
              int xpix = (xmcu * hmax + hblk) * q + x * q / 8;
              int ypix = (ymcu * vmax + vblk) * q + y * q / 8;
              int dest = (ypix * aWidth + xpix) * 4;
              int source = BLK*64*4 + (8*y+x)*4;
              if (xpix < aWidth && ypix < aHeight) {
                for (int j = 0; j < 4; j++) {
                  pixels[dest + j] = padded[source + j];
                }
              }
            }
          }
        }
      }
    }
  }
  std::cerr << "done.\n";
  
  delete this->padded;
#undef q
#undef BLK
}
