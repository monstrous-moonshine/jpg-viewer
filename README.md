# jpg-viewer

This project is my attempt to understand lossy image compression, in particular the JPEG file format. JPEG employs
a host of clever algorithms to compress image files. In particular, it

- uses a color space transformation (RGB to YCbCr) to split the image into luminance (brightness) and color
- downsamples the color components since human eyes are less sensitive to variation in color than to variation in brightness
- splits the image into 8x8 blocks, and applies the discrete cosine transform (DCT) to each block
- discards the high frequency components in each block from the frequency space representation obtained at the previous step
- lays out the 64 components in each 8x8 block in a line following a zig-zag pattern, and applies lossless compression
(Huffman encoding)

You can learn more about JPEG from [Wikipedia](https://en.wikipedia.org/wiki/JPEG) or the [official specification]
(https://www.w3.org/Graphics/JPEG/itu-t81.pdf).

## How to use
To use the program, you need to first compile it. Assuming you have cloned the repo and `cd`'d into the project
directory, compiling is as simple as typing `make` into the shell. You'll need to have [SFML](https://www.sfml-dev.org/)
installed.

Once you have compiled the program, you can invoke it this way:  
`./app [jpeg file]`

You can start with the file `Sego_lily_cm.jpg` included in the repo.

![Sego lily](Sego_lily_cm.jpg)
