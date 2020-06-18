#ifndef HUFFMAN_H
#define HUFFMAN_H

#define TABLESIZE 512
#define NONTERMINAL -1
#define INVALID 0

struct huffnode {
  int value;
  int child[2];
};

huffnode* make_hufftable(int bits[], int huffval[]);

#endif
