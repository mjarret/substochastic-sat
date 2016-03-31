/** @file  bitstring.c
 * @brief Source file for a weighted bitstring.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 3/31/16.
 */

#include <stdio.h>
#include <string.h>
#include "bitstring.h"

/**
 * Since \a blen is in units of words, this routine can be called to do the conversion automatically.
 * @param num_bits the size of bitstrings in bits.
 * @return Zero if successfull, error codes(s) if failed.
 */
int setBitLength(int num_bits){
  if ( num_bits <= 0 ) {
    return OUT_OF_BOUNDS_ERROR;
  }
  
  nbts = num_bits;
  blen = (num_bits-1)/BITS_PER_WORD + 1;
  bsize = blen*sizeof(word_t);
  return 0;
}

// Added by SPJ 3/17/16
/**
 * This prints the bitstring in the format demanded by the MaxSAT competition:
 * namely bits are numbered 1,2,3,4,... and we output -x if bit x is 0 and output x if bit x is 1.
 * @param fp is the file stream to put the output (typically stdout).
 * @param bst is the instance to be printed.
 * @return None.
 */
void printBits(FILE *fp, Bitstring bst) {
  int i;
  int val;
  
  fprintf(fp,"o %i\n", (int) bst->potential);
  fprintf(fp,"v ");
  
  for(i = 0; i < nbts; i++) {
    val = bst->node[i/BITS_PER_WORD] >> (i % BITS_PER_WORD);
    if(val%2) fprintf(fp,"%i ", i+1);
    else fprintf(fp,"%i ", -(i+1));
  }
  printf("\n");
}

/**
 * A bitstring is initialized and an array of the appropriate length is created.
 * The length will be according to the global variable \a blen.
 * @param bst_ptr points to the bitstring to be created.
 * @return Zero if successful, error code(s) if failed.
 */
int initBitstring(Bitstring *bst_ptr){
  Bitstring bst;
  
  // Create memory for the bitstring, or return failure.
  if ( (bst = (Bitstring) malloc(sizeof(struct bitstring_st))) == NULL ) {
    return MEMORY_ERROR;
  }
  
  // Create an array for storing the bits, or return failure.
  if ( (bst->node = (word_t *) calloc(blen,sizeof(word_t))) == NULL ) {
    free(bst);
    return MEMORY_ERROR;
  }
  
  // Initialize the potential to be negative, the species to be zero, and pass the created bitstring to output.
  bst->potential = -1.0;
  bst->species = 0;
  (*bst_ptr) = bst;
  
  return 0;
}

/**
 * A bitstring is destroyed and it's memory deallocated.
 * @param bst_ptr points to the instance to be destroyed.
 * @return None.
 */
void freeBitstring(Bitstring *bst_ptr){
  if ( (*bst_ptr) != NULL ) {
    if ( (*bst_ptr)->node != NULL ){
      free((*bst_ptr)->node);
    }
    free(*bst_ptr);
    (*bst_ptr) = NULL;
  }
}

/**
 * A bitstring is randomized by calls to the native LCG to C (for speed).
 * It is assumed that the randomizer has already been seeded appropriately.
 * @param bst is the instance to be randomized.
 * @return Zero.
 */
int randomBitstring(Bitstring bst){
  int i, j;
  int s = (BYTES_PER_WORD-1)/4 + 1;
  
  for (i=0; i<blen; ++i) {
    
    // Zero out the i-th word of the bitstring.
    bst->node[i] = 0ull;
    
    // Load up the i-th word of the bitstring with random bits.
    // Recall that lrand48() returns 32 bits (4 bytes) worth of random.
    for (j=0; j<s; ++j) {
      bst->node[i] |= ((word_t) lrand48()) << (32*j);
    }
  }
  
  return 0;
}

/**
 * One bitstring is copied onto another.
 * It is assumed that the memory of the target bitstring is allocated.
 * The weight is also copied.
 * @param bst_out is the location of new copy.
 * @param bst_int is the instance to be copied.
 * @return None.
 */
void copyBitstring(Bitstring bst_out, Bitstring bst_in){
  
  // Copy each element of the array.
//  for (i=0; i<blen; ++i)
//    bst_out->node[i] = bst_in->node[i];
  memcpy(bst_out->node, bst_in->node, bsize);
  
  // Copy the weight and species tag.
  bst_out->potential = bst_in->potential;
  bst_out->species = bst_in->species;
}

/**
 * A random bit is flipped in the input bitstring, and result is copied to the output.
 * It is assumed that the memory of the target bitstring is allocated.
 * The weight is also copied.
 * The index (1-up) of the bit flipped is returned, with a sign: if negative the flip went 1->0, if positive 0->1.
 * @param bst_out is the location of new bitstring.
 * @param bst_int is the instance where a bit is flipped.
 * @return Signed index of bit that was flipped.
 */
int randomBitFlip(Bitstring bst_out, Bitstring bst_in){
  int i;
  word_t j;
  
  // Copy bitstring.
  copyBitstring(bst_out, bst_in);
  
  // Choose a random bit.
  i = lrand48() % nbts;
  
  // Find the value of the bit.
  j = (bst_out->node[i/BITS_PER_WORD] >> (i % BITS_PER_WORD)) & 1;
  
  // Flip that bit.
  bst_out->node[i/BITS_PER_WORD] ^= ((word_t) 1) << (i % BITS_PER_WORD);
  
  // Return the bit that was flipped.
  return (1-2*((int) j))*(i+1);
}


