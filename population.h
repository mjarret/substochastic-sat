/** @file  population.h
 * @brief Header file for a population.
 *
 * Created by Brad Lackey on 3/30/16. Last modified 3/31/16.
 */

#ifndef population_h
#define population_h

#include <stdio.h>
#include <stdlib.h>

#include "macros.h"
#include "bitstring.h"
#include "sat.h"

int arraysize; ///< The integer \a arraysize is the amount of memory assigned to the population.
int nspecies; /// The integer \a nspecies is the number of subpopulations in the population.


/// The underlying type for a population.
/**
 * This data type also holds a reference to the SAT instance.
 */
struct population_st;
typedef struct population_st * Population;

struct population_st {
  SAT sat;             ///< Reference to the underlying SAT problem.
  DSAT ds;             ///< Reference to the SAT derivative for updates.
  int size;            ///< The overall current population size.
  int *psize;          ///< The current size of each subpopulation.
  Bitstring *walker;   ///< Array of bitstrings that form the population.
  Bitstring winner;    ///< A copy of the best walker yet found.
  double *avg_v;       ///< Array of average potentials for each subpopulation.
  double *max_v;       ///< Maximum potential for each subpopulation.
  double *min_v;       ///< Minimum potential for each subpopulation.
};

// Memory management routines.
int initPopulation(Population *Pptr, SAT sat);   ///< Create a population from its SAT instance.
void freePopulation(Population *Pptr);           ///< Deallocation routine for a population.
int randomPopulation(Population P, int size);    ///< Initialize a population of size \a size with random walkers.

#endif /* population_h */
