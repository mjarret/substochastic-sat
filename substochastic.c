/** @file  substochastic.c
 * @brief Main file for Substochastic Monte Carlo.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 5/18/16.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "macros.h"
#include "bitstring.h"
#include "sat.h"
#include "population.h"

//On machines with very old versions of glibc (e.g. the Raritan cluster)
//we need to define gnu_source in order to avoid warnings about implicit
//declaration of getline() and round().
#define _GNU_SOURCE

extern int blen;
extern int nbts;
extern int arraysize;
extern int problem_type;
extern potential_t topweight;

static int popsize,runmode;
static double weight, end_weight, runtime, runstep;
static potential_t optimal;




void update(double a, double b, double mean, Population P, int parity);
int parseCommand(int argc, char **argv, Population *Pptr);


int main(int argc, char **argv){
  int parity, try, err;
  double mean;
  double a, b, t, dt;
  Population pop;
  potential_t local_min, min = -1;      //the best minimum from different trials
  Bitstring solution;   //the corresponding bitstring
  clock_t beg, end;     //for code timing
  double time_spent;    //for code timing
  
#if TRACK_GLOBAL_BIASES
  int i;
  word_t u;
#endif
  
  
  beg = clock();
  
  if ( (err = parseCommand(argc, argv, &pop)) ){
    return err;
  }
  
  end = clock();
  time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
  printf("c Problem loaded: %f seconds\n", time_spent);

  if ( (err = initBitstring(&solution)) ){
    fprintf(stderr, "Could not initialize answerspace.\n");
    return err;
  }
  
  
  randomPopulation(pop,popsize);
  end = clock();
  if ( pop->winner->potential < topweight ) {
    printf("o %ld\n", pop->winner->potential);
    fflush(stdout);
    printBits(stdout, pop->winner);
    fflush(stdout);
    time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
    printf("c Walltime: %f seconds, 0 loops\n", time_spent);
  }
  fflush(stdout);
  min = pop->winner->potential;
  if (min <= optimal) exit(0);
  
  try = 1;
  while (1) {
    
    t = 0.0;
    parity = 0;
    randomPopulation(pop,popsize);
    local_min = min;

    while (t < runtime) {
      
      // The annealing schedule
      a = weight*(1.0 - t/runtime); // Turned weight into percent -- Michael 3/30/16
      b = (t/runtime);
      
      mean = pop->avg_v + (pop->max_v - pop->min_v)*(popsize - pop->psize)/(2.0*popsize);
      
      if ( (pop->max_v - mean) > (mean - pop->min_v) )
        dt = 0.9/(a + b*(pop->max_v - mean));
      else
        dt = 0.9/(a + b*(mean - pop->min_v));
      if (t + dt > runtime)
        dt = runtime - t;
      
      update(a*dt, b*dt, mean, pop, parity);
      
      end = clock();
      time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
      
      if ( pop->winner->potential < local_min ) {
        local_min = pop->winner->potential;
        if ( local_min < topweight ) {
          printf("o %ld\n", local_min);
          printf("c Walltime: %f seconds, %d loop(s)\n", time_spent, try);
          fflush(stdout);
        }
        if (local_min <= optimal) {
          break;
        }
      }
      if ( time_spent > 290 )
        break;

      t += dt;
      parity ^= 1;
      
    }
    
    if ( local_min < min ) {
      min = local_min;
      copyBitstring(solution, pop->winner);
      if ( min < topweight ) {
        printBits(stdout, solution);
        fflush(stdout);
      }
      if (min <= optimal) {
        sleep(1);
        return 0;
      }
    }
    
    if ( time_spent > 290 ){
      return 1;
    }
    
    if ( time_spent > 120 ){
      
      weight = pop->sat->total_weight/150000.0;
      runtime = 8000000.0;
      
    } else {
      
      weight = 0.95*weight + 0.05*end_weight;
      runtime += runstep;
      
    }
    

#if TRACK_GLOBAL_BIASES
    
    for (i=0; i<pop->sat->num_vars; ++i) {
      pop->sat->global_bias[i] *= REMIX_PERCENTAGE;
      u = pop->winner->node[i/VARIABLE_NUMB_BITS] >> (i % VARIABLE_NUMB_BITS);
      pop->sat->global_bias[i] += (1-REMIX_PERCENTAGE)*(UPDATE_RELAXATION - (2*UPDATE_RELAXATION -1)*(u%2));
    }
    
#endif
    
    ++try;
    randomPopulation(pop,popsize);
  }
  
  freeBitstring(&solution);

  return 0;
}


int parseCommand(int argc, char **argv, Population *Pptr){
  SAT sat;
  int seed;
  FILE *fp;
  Population pop;
  
  if ( argc < 2 || argc > 4 ) {
    fprintf(stderr, "Usage: %s <instance.cnf> \n",argv[0]);
    fprintf(stderr, "Usage: %s <instance.cnf> [<target optimum> [<seed>]]\n",argv[0]);
    return 2;
  }
  
  printf("c ------------------------------------------------------\n");
  printf("c ----------------------------------------------------------\n");
  printf("c Substochastic Monte Carlo, version 1.0, 2016\n");
  printf("c Michael Jarret, Stephen Jordan, and Brad Lackey\n");
  printf("c Joint Center for Quantum Information and Computer Science\n");
  printf("c University of Maryland, College Park.\n");
  printf("c ----------------------------------------------------------\n");
  printf("c Input: %s\n", argv[1]);
  
  if ( (fp = fopen(argv[1], "r")) == NULL ){
    fprintf(stderr,"Could not open file %s\n",argv[1]);
    return IO_ERROR;
  }
  
  if ( loadDIMACSFile(fp,&sat) ){
    fprintf(stderr,"Error reading in DIMACS SAT file %s\n",argv[1]);
    return IO_ERROR;
  }
  
  fclose(fp);
  
  setBitLength(sat->num_vars);
  
  printf("c Bits: %d\n", nbts);
  printf("c Clauses (after tautology removal): %d\n", sat->num_clauses);
  printf("c Problem type: %d\n", problem_type);
  
  if ( problem_type == UNKNOWN ) {
    weight = sat->total_weight/5000.0;
    end_weight = 0.01;
    runtime = sat->num_vars*sat->num_vars/100;
    runstep = 2*runtime;
    popsize = 1024;
    if ( NUM_VARIABLE_WORDS*vlen*tlen*sizeof(word_t) + NUM_CLAUSE_WORDS*clen*sizeof(int) < (1<<30) )
      runmode = 2;
    else
      runmode = 1;
  }
  
  if ( problem_type == UNWEIGHTED_2_SAT ){
    
    // Parameters finalized: 13 April 2016 (bclackey)
    weight = sat->total_weight/5000.0;
    end_weight = sat->total_weight/5500.0;
    runtime = exp(0.022*sat->num_vars + 5.9);
    runstep = exp(0.022*sat->num_vars + 4.6);
    popsize = 16;
    runmode = 1;
    
  }
  
  if (  problem_type == PARTIAL_2_SAT ){
    
    // Parameters finalized: 14 April 2016 (bclackey)
    weight = sat->total_weight/60000.0;
    end_weight = sat->total_weight/100000.0;
    runtime = exp(0.018*sat->num_vars + 3.7);
    runstep = exp(0.012*sat->num_vars + 4.3);
    popsize = 16;
    runmode = 1;
    
  }
  
  if ( problem_type == WEIGHTED_2_SAT ){
    
    // Parameters finalized: 13 April 2016 (bclackey)
    weight = sat->total_weight/5000.0;
    end_weight = sat->total_weight/5500.0;
    runtime = exp(0.022*sat->num_vars + 5.9);
    runstep = exp(0.022*sat->num_vars + 4.6);
    popsize = 16;
    runmode = 1;
    
    
  }
  
  if ( problem_type == UNWEIGHTED_3_SAT ) {
    
    // Parameters finalized: 13 April 2016 (bclackey)
    weight = sat->total_weight/5000.0;
    end_weight = sat->total_weight/10000.0;
    runtime = exp(0.035*sat->num_vars + 6.1);
    runstep = exp(0.030*sat->num_vars + 4.4);
    popsize = 16;
    runmode = 1;
    
  }
  
  if ( problem_type == PARTIAL_3_SAT ) {
    
    // Parameters finalized: 14 April 2016 (bclackey)
    weight = sat->total_weight/16000.0;
    end_weight = sat->total_weight/200000.0;
    runtime = exp(0.031*sat->num_vars + 4.4);
    runstep = exp(0.025*sat->num_vars + 3.9);
    popsize = 16;
    runmode = 1;
    
  }
  
  if ( problem_type == WEIGHTED_3_SAT ){
    
    // Parameters finalized: 13 April 2016 (bclackey)
    weight = sat->total_weight/5000.0;
    end_weight = sat->total_weight/10000.0;
    runtime = exp(0.028*sat->num_vars + 6.2);
    runstep = exp(0.028*sat->num_vars + 4.5);
    popsize = 16;
    runmode = 1;
    
  }
  
  if ( (problem_type == UNWEIGHTED_4_SAT) || (problem_type == WEIGHTED_4_SAT) ) {
    
    // Parameters last modified: 13 April 2016 (bclackey)
    weight = sat->total_weight/8000.0;
    end_weight = 0.01;
    runtime = exp(0.032*sat->num_vars + 9.3);
    runstep = exp(0.036*sat->num_vars + 7.9);
    popsize = 16;
    runmode = 2;
  }
  
  
  if ( argc >= 3 ) {
    
    optimal = atoi(argv[2]);
    
    if ( argc == 4 )
      seed = atoi(argv[3]);
    else
      seed = time(0);

  } else {
      
    optimal = 0;
    seed = time(0);
    
  }

  
  // Break out of attempt to put way too large a problem into the system.
  if ( runtime > 8000000 ) {
    runtime = 8000000;
    runstep = 0;
    weight = 0.1;
    end_weight = 0.01;
  }
  if ( runtime < 500 ){
    runtime = 500;
    runstep = 100;
  }
  
  //  printf("c Population size: %d\n", popsize);
  //  printf("c Starting runtime: %.0f\n", runtime);
  //  printf("c Runtime step per loop: %.0f\n", runstep);
  //  printf("c Step weight: %.3f\n", weight);
  printf("c Target potential: %ld\n", optimal);
  printf("c Top potential: %ld\n", topweight);
  arraysize = 2000;
  
  srand48(seed);
  printf("c Seed: %i\n", seed);

  if ( initPopulation(&pop, sat, runmode) ) {
    fprintf(stderr,"Could not initialize potential.\n");
    return MEMORY_ERROR;
  }
  
  *Pptr = pop;
  
  return 0;
}


void update(double a, double b, double mean, Population P, int parity){
  int i,j,k;
  double p,e;
  potential_t min = P->sat->num_clauses;
  double avg = 0.0;
  potential_t max = -(P->sat->num_clauses);
  int old = parity*arraysize;
  int new = (1-parity)*arraysize;
  
  for (i=j=0; i<P->psize; ++i) {   // Loop over each walker (i) and set target position (j) to zero.
    p = drand48();
    
    // First potential event: walker steps.
    if ( p < a ) {
      k = randomBitFlip(P->walker[new+j], P->walker[old+i]);
      if ( P->ds != NULL )
        e = P->walker[old+i]->potential + ((k>0)-(k<0))*getPotential(P->walker[new+j],P->ds->der[abs(k)-1]);
      else {
        if ( P->tbl != NULL )
          e = getPotential2(P->walker[new+j],P->tbl);
        else
          e = getPotential(P->walker[new+j],P->sat);
      }

      P->walker[new+j]->potential = e;
      if ( e < min ){
        min = e;
        if ( e < P->winner->potential )
          copyBitstring(P->winner, P->walker[new+j]);
      }
      if ( e > max ) max = e;
      avg += e;
      ++j;
      continue;
    }
    p -= a;
    
    // Second potential event: walker spawns/dies.
    e = b*(mean - P->walker[old+i]->potential);
    if ( p < e ) {  // Spawn.
      copyBitstring(P->walker[new+j], P->walker[old+i]);
      copyBitstring(P->walker[new+j+1], P->walker[old+i]);
      e = P->walker[old+i]->potential;
      if ( e < min ) min = e;
      if ( e > max ) max = e;
      avg += 2*e;
      j += 2;
      continue;
    }
    
    if ( p < -e ) { // Die.
      continue;
    }
    
    // Third potential event: walker stays.
    copyBitstring(P->walker[new+j], P->walker[old+i]);
    e = P->walker[old+i]->potential;
    if ( e < min ) min = e;
    if ( e > max ) max = e;
    avg += e;
    ++j;
  }
  
  P->psize = j;
  P->avg_v = avg/j;
  P->min_v = min;
  P->max_v = max;
}
