#ifndef JPEG_H
#define JPEG_H

#include "huffman.h"

#define SOF 0xc0
#define DHT 0xc4
#define DAC 0xcc
#define RST 0xd0
#define SOI 0xd8
#define EOI 0xd9
#define SOS 0xda
#define DQT 0xdb
#define DRI 0xdd
#define APP 0xe0
#define JPG 0xf0
#define COM 0xfe

class image_t {
  // private variables
  unsigned char* buffer;
  int fileSize;
  int pos = 0;
  
  // display dimensions
  int dispWidth, dispHeight;
  int quality;
  
  // image parameters
  int nmcus, size_mcu = 0;
  int hmax = 0, vmax = 0;
  int Ci[3];
  int Hi[3], Vi[3];
  int Tq[3];
  int rst = 0;
  
  // intermediate buffer
  int* data;
  unsigned char* padded;
  
  // tables
  huffnode* htables[4];
  int* qtables[2];
  
  // private functions
  void print_markers();
  void read_jpeg();
  void read_marker(int marker_type);
  void read_tables();
  void read_hufftable(int);
  void read_quantable(int);
  void read_rst();
  void read_frame();
  void read_scan();
  void dequantize();
  void unzigzag();
  void idct();
  void color_convert();
  void pixel_rearrange();
  
public:
  int height, width;
  int aHeight, aWidth; // adjusted height and width
  unsigned char* pixels;
  image_t(const char* filePath, int width, int height);
  ~image_t();
};

#endif
