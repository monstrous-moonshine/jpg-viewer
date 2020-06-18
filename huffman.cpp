#include <iostream>
#include <memory>
#include "huffman.h"

void traverse(huffnode* hufftable, int idx, int level) {
    huffnode* entry = &hufftable[idx];
    for (int i = 0; i < level; i++) std::cerr << " ";
    std::cerr << "[" << idx << "] " << entry->value;
    std::cerr << " [" << entry->child[0] << " " << entry->child[1] << "]\n";
    if (entry->child[0] != INVALID) {
        traverse(hufftable, entry->child[0], level + 1);
    }
    if (entry->child[1] != INVALID) {
        traverse(hufftable, entry->child[1], level + 1);
    }
}

// we pad bits with a byte in front so the indices
// correspond to number of bits
huffnode* make_hufftable(int bits[], int huffval[]) {
    int total = 0;
    
    // tally up bits
    for (int nbits = 1; nbits <= 16; nbits++) {
        total += bits[nbits];
    }
    
    // generate huffcode
    std::unique_ptr<int[]> huffcode(new int[total]);
    for (int nbits = 1, k = 0, code = 0; nbits <= 16; nbits++) {
        for (int i = 0; i < bits[nbits]; i++) {
            huffcode[k++] = code++;
        }
        code <<= 1;
    }
    
    // initialize hufftable
    huffnode* hufftable(new huffnode[TABLESIZE]);
    for (int i = 0; i < TABLESIZE; i++) {
        hufftable[i].value = NONTERMINAL;
        hufftable[i].child[0] = INVALID;
        hufftable[i].child[1] = INVALID;
    }
  
    // generate hufftable
    for (int nbits = 1, k = 0, next = 1; nbits <= 16; nbits++) {
        for (int i = 0; i < bits[nbits]; i++) {
        for (int pos = 0, idx = 0; pos < nbits; pos++) {
            int c = (huffcode[k] >> (nbits - 1 - pos)) & 1;
            if (hufftable[idx].child[c] == INVALID) {
                hufftable[idx].child[c] = next;
                idx = next;
                next++;
            } else {
                idx = hufftable[idx].child[c];
            }
            
            // guard against table overflow
            if (next >= TABLESIZE) {
                std::cerr << "Table limit exceeded.\n";
                exit(EXIT_FAILURE);
            }
        }
        
        // assign value
        hufftable[next - 1].value = huffval[k++];
        }
    }
    // traverse(hufftable, 0, 0);
    return hufftable;
}
