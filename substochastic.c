/** @file  substochastic.c
 * @brief Main file for Substochastic Monte Carlo.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 3/30/16.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "macros.h"
#include "bitstring.h"
#include "sat.h"
#include "population.h"

//On machines with very old versions of glibc (e.g. the Raritan cluster)
//we need to define gnu_source in order to avoid warnings about implicit
//declaration of getline() and round().
#define _GNU_SOURCE

#define TEST 1

extern int nbts;
extern int arraysize;
extern int nspecies;
extern size_t bsize;

static int popsize;

static double e; // Temporary register


void update(double a, double b, double *mean, Population P, int parity);
int parseCommand(int argc, char **argv, Population *Pptr, double *weight, double *runtime, int *trials, double *starttime);


/*Added by SPJ 3/17
 *Print optimal bitstring out of the current population. If there are multiple
 *optima with the same potential, just print one of them. Returns the index of
 *the optimum printed.
 */
/*int printOpt(Population P) {
  int i;
  double potential;
  i = -1;
  do {
    i++;
    potential = getPotential(P->walker[i], P->sat);
    if(potential == P->min_v) printBits(P->walker[i]);
  }while(potential != P->min_v);
  
  return i;
}*/

/*Added by SPJ 3/17
 *Returns the average number of variables per clause.
 */
double avgLength(SAT instance) {
  int c;
  int total;
  total = 0;
  for(c = 0; c < instance->num_clauses; c++) total += instance->clause_length[c];
  return (double)total/(double)instance->num_clauses;
}

/*Added by SPJ 3/17
 *Automatically chooses parameters, such as runtime and number of trials.
 */
void autoparam(double varsperclause, int vars, double *weight, double *runtime, int *trials, double *starttime) {
  int k;
  k = (int)round(varsperclause);
  //*weight = 1.0;
  *weight = 100.0; // Changed since weight now percent -- Michael 3/30/16
  *trials = 5; //default
  *starttime = 0.0; // default -- Michael 3/30/16
  if(k == 2) {
    //These values are determined using the experimental results in
    //summary_2sat120v.txt and summary_2sat200v.txt.
    *trials = 10;
    *runtime = 3.2E4;
  }
  if(k == 3) {
    //These values are determined using the experimental results in
    //summary_3sat70v.txt and summary3sat110vlong.txt.
    *trials = 5;
    *runtime = 223.0*exp(0.07*vars);
  }
  if(k == 4 || vars > 200) {
    //This overrides the above.
    //Just go for broke and use up all the time.
    *trials = 5;
    *runtime = 3E6;
  }
}

int main(int argc, char **argv){
  int i, parity, trials, try, err;
  double weight, runtime, starttime, max_r;
  double *mean;
  double a, b, t, dt;
  Population pop;
  double min = -1.0;    //the best minimum from different trials
  double local_min = -1.0;
  Bitstring solution;   //the corresponding bitstrings
  clock_t beg, end;     //for code timing
  double time_spent;    //for code timing
  
  beg = clock();

  if ( (err = parseCommand(argc, argv, &pop, &weight, &runtime, &trials, &starttime)) ){
    return err;
  }
  
  if ( (mean = (double *) malloc(nspecies*sizeof(double))) == NULL ) {
    fprintf(stderr, "Could not initialize means.\n");
    return MEMORY_ERROR;
  }
  
  if ( (err = initBitstring(&solution)) ){
    fprintf(stderr, "Could not initialize answerspace.\n");
    return err;
  }
  
  for(try = 0; try < trials; try++) {
    t = starttime;
    parity = 0;
    randomPopulation(pop,popsize);
    
    while (t < runtime) {
      a = weight*(1.0 - t/runtime)/100.0; // Turned weight into percent -- Michael 3/30/16
      b = (t/runtime);
      
      max_r = 0.0;
      for (i=0; i<nspecies; ++i){
        mean[i] = pop->avg_v[i] + (pop->max_v[i] - pop->min_v[i])*(popsize - pop->psize[i])/(2*popsize);
        if ( (pop->max_v[i] - mean[i]) > (mean[i] - pop->min_v[i]) ) {
          if ( (pop->max_v[i] - mean[i]) > max_r )
            max_r = pop->max_v[i] - mean[i];
        } else {
          if ( (mean[i] - pop->min_v[i]) > max_r )
            max_r = mean[i] - pop->min_v[i];
        }
      }
      
      dt = 0.9/(a + b*max_r);
      
      
/*      if ( (local_min<0) || (pop->winner->potential < local_min)) {
        local_min = pop->winner->potential;
        printf("%f: ",t/runtime); for (i=0; i<nspecies; ++i) printf("(%4.0f, %4.1f, %4.0f) ",pop->max_v[i],mean[i],pop->min_v[i]); printf("\n");
      }
*/
      
      if (t + dt > runtime)
        dt = runtime - t;
      
      update(a*dt, b*dt, mean, pop, parity);
      
      t += dt;
      parity ^= 1;
    }
    
#if TEST
    printBits(stdout, pop->winner);
    if ((min<0) || (pop->winner->potential < min)) {
      min = pop->winner->potential;
      copyBitstring(solution, pop->winner);
    }
#else
    end = clock();
    time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
    if ((min<0) || (pop->winner->potential < min)) {
      min = pop->winner->potential;
      copyBitstring(solution, pop->winner);
      printBits(stdout, pop->winner);
      printf("c Walltime: %f seconds\n", time_spent);
      fflush(stdout);
    } else {
      if ( time_spent > 240 ) {
        break;
      }
    }
    
#endif
  }

#if 1
//  printf("c Final answer: \n");
//  printBits(stdout, solution);
  freeBitstring(&solution);
  end = clock();
  time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
  printf("c Walltime: %f seconds\n", time_spent);
#endif
  return 0;
}


void update(double a, double b, double *mean, Population P, int parity){
  int i,j,k,l;
  double p;
  Bitstring in, out;
  int old = parity*arraysize;
  int new = (1-parity)*arraysize;
  
  for (i=0; i<nspecies; ++i) {
    P->psize[i] = 0;
    P->avg_v[i] = 0.0;
    P->max_v[i] = -(P->sat->num_clauses);
    P->min_v[i] = P->sat->num_clauses;
  }
  
  
  for (i=j=0; i<P->size; ++i) {   // Loop over each walker (i) and set target position (j) to zero.
    p = drand48();
    in = P->walker[old+i];
    out = P->walker[new+j];
    l = P->walker[old+i]->species; // This is the species where the walker will go!
    
    
    // First potential event: walker steps.
    if ( p < a ) {
      k = randomBitFlip(out, in);
      e = in->potential + ((k>0)-(k<0))*getPotential(out,P->ds->der[abs(k)-1]);
      out->potential = e;
      if ( e < P->min_v[l] ){
        P->min_v[l] = e;
        if ( e < P->winner->potential )
          copyBitstring(P->winner, out);
      }
      if ( e > P->max_v[l] ) P->max_v[l] = e;
      P->avg_v[l] += e;
      ++(P->psize[l]);
      ++j;
      continue;
    }
    p -= a;
   
    e = b*(mean[in->species] - in->potential);
    // Second potential event: walker spawns/dies.
    if ( p < e ) {
//      copyBitstring(out, in);
      memcpy(out->node, in->node, bsize);
      e = out->potential = in->potential;
      out->species = in->species;
//      copyBitstring(P->walker[new+j+1], in);
      out = P->walker[new+j+1];
      memcpy(out->node, in->node, bsize);
      out->potential = in->potential;
      if ( e < P->min_v[l] ) P->min_v[l] = e;
      if ( e > P->max_v[l] ) P->max_v[l] = e;
      P->avg_v[l] += 2*e;
      P->psize[l] += 2;
      j += 2;
      continue;
      }
    
    if ( p < -e ) { // Case of dying.
      continue;
    }
    
    // Third potential event: walker stays.
//    copyBitstring(out, in);
    memcpy(out->node, in->node, bsize);
    e = out->potential = in->potential;
    out->species = in->species;
    if ( e < P->min_v[l] ) P->min_v[l] = e;
    if ( e > P->max_v[l] ) P->max_v[l] = e;
    P->avg_v[l] += e;
    ++(P->psize[l]);
    ++j;
  }
  
  P->size = j;
  
  for (i=0; i<nspecies; ++i) P->avg_v[i] /= P->psize[i];
}



int parseCommand(int argc, char **argv, Population *Pptr, double *weight, double *runtime, int *trials, double *starttime){
  SAT sat;
  int seed;
  double varsperclause; //average number of variables per clause
  FILE *fp;
  Population pop;
  
  if (argc != 8 && argc != 2) {
    fprintf(stderr, "Usage: %s instance.cnf [<step weight> <runtime> <species size> <number species> <trials> <start time>]\n",argv[0]);
    return 2;
  }
  
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
  varsperclause = avgLength(sat);

  if ( argc == 2 )
    autoparam(varsperclause, nbts, weight, runtime, trials, starttime);
  
  if( argc == 8 ) {
    *weight = (double) atoi(argv[2]);
    sscanf(argv[3], "%lf", runtime); //allows scientific notation unlike atoi
    popsize = atoi(argv[4]);
    nspecies = atoi(argv[5]);
    *trials = atoi(argv[6]);
    *starttime = atoi(argv[7]);
  } else {
    popsize = 128;
    nspecies = 1;
  }
  
  arraysize = 20*nspecies*popsize;
  
  if ( initPopulation(&pop, sat) ) {
    fprintf(stderr,"Could not initialize potential.\n");
    return MEMORY_ERROR;
  }

  seed = time(0);
//  seed = 0; // for testing
  srand48(seed);
  printf("c ------------------------------------------------------\n");
  printf("c Substochastic Monte Carlo, version 1.0                \n");
  printf("c Brad Lackey, Stephen Jordan, and Michael Jarret, 2016.\n");
  printf("c ------------------------------------------------------\n");
  printf("c Input: %s\n", argv[1]);
  printf("c Bits: %d\n", nbts);
  printf("c Clauses (after tautology removal): %d\n", pop->sat->num_clauses);
  printf("c Step weight: %f\n", *weight);
  printf("c Population size: %d (=%dx%d)\n", nspecies*popsize, popsize, nspecies);
  printf("c Runtime: %e\n", *runtime);
  printf("c Start time: %e\n", *starttime);
  printf("c Trials: %i\n", *trials);
  printf("c Variables per clause: %f\n", varsperclause);
  printf("c Seed: %i\n", seed);

  *Pptr = pop;
  
  return 0;
}




