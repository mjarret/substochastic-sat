/** @file  population.c
 * @brief Source file for a population.
 *
 * Created by Brad Lackey on 3/30/16. Last modified 3/31/16.
 */

#include <stdio.h>
#include <stdlib.h>
#include "population.h"



/**
 * This creates a population from the passed SAT instance, and computes and stores the derivative of that SAT instance too.
 * The amount of memory assigned to the population is an global integer \a arraysize.
 * Note the population size is dynamics, and so enough memory must be assigned.
 * WARNING: there are currently no checks to ensure the population does not die off entirely, or overflow the array!
 * However in the current implementation of the \a update routine, it is highly unlikely that the population size grows/shrinks by more that a few percent of its initial size.
 * @param Pptr points to the population to be created.
 * @param sat is the underlying SAT instance.
 * @return Zero if successful, error code(s) if failed.
 */
int initPopulation(Population *Pptr, SAT sat){
  int i,err;
  Population P;
  
  if ( (P = (Population) malloc(sizeof(struct population_st))) == NULL) {
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  P->sat = sat;
  createSATDerivative(&(P->ds),sat);
  
  if ( (P->walker = (Bitstring *) malloc((2*arraysize)*sizeof(Bitstring))) == NULL ){
    freePopulation(&P);
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  if ( (err = initBitstring(&(P->winner))) ) {
    freePopulation(&P);
    *Pptr = NULL;
    return err;
  }
  
  for (i=0; i<(2*arraysize); ++i) {
    if ( (err = initBitstring(P->walker + i)) ) {
      freePopulation(&P);
      *Pptr = NULL;
      return err;
    }
  }
  
  if ( (P->psize = (int *) calloc(nspecies,sizeof(int))) == NULL ) {
    freePopulation(&P);
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  if ( (P->avg_v = (double *) calloc(nspecies,sizeof(double))) == NULL ) {
    freePopulation(&P);
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  if ( (P->max_v = (double *) calloc(nspecies,sizeof(double))) == NULL ) {
    freePopulation(&P);
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  if ( (P->min_v = (double *) calloc(nspecies,sizeof(double))) == NULL ) {
    freePopulation(&P);
    *Pptr = NULL;
    return MEMORY_ERROR;
  }
  
  *Pptr = P;
  return 0;
}


/**
 * @param Pptr points to the population instance to be destroyed.
 * @return None.
 */
void freePopulation(Population *Pptr){
  int j;
  if ( *Pptr != NULL ) {
    if ( (*Pptr)->sat != NULL )
      freeSAT(&((*Pptr)->sat));
    if ( (*Pptr)->ds != NULL )
      freeSATDerivative(&((*Pptr)->ds));
    if ( (*Pptr)->winner != NULL )
      freeBitstring(&((*Pptr)->winner));
    if ( (*Pptr)->walker != NULL ) {
      for (j=0; j<(2*arraysize); ++j)
        freeBitstring((*Pptr)->walker + j);
      free((*Pptr)->walker);
    }
    if ( (*Pptr)->psize != NULL )
      free((*Pptr)->psize);
    if ( (*Pptr)->avg_v != NULL )
      free((*Pptr)->avg_v);
    if ( (*Pptr)->max_v != NULL )
      free((*Pptr)->max_v);
    if ( (*Pptr)->min_v != NULL )
      free((*Pptr)->min_v);
    *Pptr = NULL;
  }
}


/**
 * This routine initializes an already created population with random walkers.
 * The global parameter \a nspecies must be set before calling this routine.
 * The potential of each walker is computed and stored, so updates can be performed using derivatives.
 * The walker with the best potential is stored off as the winner.
 * Note: this initial distribution is uniform, which is the ground state for the hypercube Laplacian.
 * If one changes the driving Hamiltonian, then this routine will likely need to be modified.
 * @param P is the population to be initialized.
 * @param size is the initial size of each subpopulation.
 * @return Zero.
 */
int randomPopulation(Population P, int size){
  int i, j, argmin;
  double e, avg, min, max, global_min;
  
  randomBitstring(P->walker[0]);
  P->walker[0]->species = 0;
  e = getPotential(P->walker[0], P->sat);
  P->walker[0]->potential = e;
  avg = e;
  min = e;
  global_min = e;
  argmin = 0;
  max = e;
  
  for (j=1; j<size; ++j){
    randomBitstring(P->walker[j]);
    P->walker[j]->species = 0;
    e = getPotential(P->walker[j], P->sat);
    P->walker[j]->potential = e;
    avg += e;
    if ( e < min ){
      min = e;
      global_min = e;
      argmin = j;
    }
    if ( e > max )
      max = e;
  }
  P->psize[0] = size;
  P->avg_v[0] = avg/size;
  P->max_v[0] = max;
  P->min_v[0] = min;
  
  
  for (i=1; i<nspecies; ++i) {
    randomBitstring(P->walker[i*size]);
    P->walker[i*size]->species = i;
    e = getPotential(P->walker[i*size], P->sat);
    P->walker[i*size]->potential = e;
    avg = e;
    min = e;
    max = e;
    if ( e < global_min ) {
      global_min = e;
      argmin = i*size;
    }
    
    for (j=1; j<size; ++j) {
      randomBitstring(P->walker[i*size+j]);
      P->walker[i*size+j]->species = i;
      e = getPotential(P->walker[i*size+j], P->sat);
      P->walker[i*size+j]->potential = e;
      avg += e;
      if ( e < min ){
        min = e;
        if ( e < global_min ) {
          global_min = e;
          argmin = i*size+j;
        }
      }
      if ( e > max )
        max = e;
    }
    P->psize[i] = size;
    P->avg_v[i] = avg/size;
    P->max_v[i] = max;
    P->min_v[i] = min;
  }
  
  P->size = nspecies*size;
  
  copyBitstring(P->winner, P->walker[argmin]);
  
  return 0;
}

