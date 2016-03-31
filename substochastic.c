/** @file  substochastic.c
 * @brief Main file for Substochastic Monte Carlo.
 *
 * Created by Brad Lackey on 3/14/16. Last modified 3/30/16.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "macros.h"
#include "bitstring.h"
#include "sat.h"
#include "population.h"

//On machines with very old versions of glibc (e.g. the Raritan cluster)
//we need to define gnu_source in order to avoid warnings about implicit
//declaration of getline() and round().
#define _GNU_SOURCE

extern int nbts;
extern int arraysize;
extern int nspecies;

static int popsize;


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
  int i,parity, trials, try, err;
  double weight, runtime, starttime, max_r;
  double *mean;
  double a, b, t, dt;
  Population pop;
  int *mins;            //the minima from different trials
  Bitstring *solutions; //the corresponding bitstrings
  int best_try;         //the trial in which the best minimum was found
  int initfail;         //a flag for memory allocation failure
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
  
  mins = (int *)malloc(trials*sizeof(int));
  solutions = (Bitstring *)malloc(trials*sizeof(Bitstring));
  if(mins == NULL || solutions == NULL) {
    fprintf(stderr, "Could not initialize trials.\n");
    return MEMORY_ERROR;
  }
  initfail = 0;
  for(try = 0; try < trials && !initfail; try++) initfail = initBitstring(&solutions[try]);
  if(initfail) {
    fprintf(stderr, "Could not initialize answerspace.\n");
    return MEMORY_ERROR;
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
        mean[i] = pop->avg_v[i] + (pop->max_v[i] - pop->min_v[i])*(popsize - pop->psize[i])/((pop->max_v[i] + pop->min_v[i])*popsize);
        if ( (pop->max_v[i] - mean[i]) > (mean[i] - pop->min_v[i]) ) {
          if ( (pop->max_v[i] - mean[i]) > max_r )
            max_r = pop->max_v[i] - mean[i];
        } else {
          if ( (mean[i] - pop->min_v[i]) > max_r )
            max_r = mean[i] - pop->min_v[i];
        }
      }
      
      dt = 0.9/(a + b*max_r);
//      printf("%f %f %d %d (%f %f %f)\n",t,dt,pop->size,pop->psize[0],pop->max_v[0],mean[0],pop->min_v[0]);
      if (t + dt > runtime)
        dt = runtime - t;
      
      update(a*dt, b*dt, mean, pop, parity);
      
      t += dt;
      parity ^= 1;
    }
//    printf("o %i\n",(int)pop->min_v);
//    printf("v ");
//    optindex = printOpt(pop);
    
    printBits(stdout, pop->winner);
    

    
    mins[try] = (int)pop->winner->potential;
//    copyBitstring(solutions[try], pop->walker[optindex]);
    copyBitstring(solutions[try], pop->winner);
  }
  best_try = 0;
  for(try = 1; try < trials; try++) if(mins[try] < mins[best_try]) best_try = try;
  printf("c Final answer: \n");
  printBits(stdout, solutions[best_try]);
  for(try = 0; try < trials; try++) freeBitstring(&solutions[try]);
  free(solutions);
  free(mins);
  end = clock();
  time_spent = (double)(end - beg)/CLOCKS_PER_SEC;
  printf("c Walltime: %f seconds\n", time_spent);
  return 0;
}


void update(double a, double b, double *mean, Population P, int parity){
  int i,j,k,l;
  double p,e;
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
    l = P->walker[old+i]->species; // This is the species where the walker will go!
    
    
    // First potential event: walker steps.
    if ( p < a ) {
      k = randomBitFlip(P->walker[new+j], P->walker[old+i]);
      e = P->walker[old+i]->potential + ((k>0)-(k<0))*getPotential(P->walker[new+j],P->ds->der[abs(k)-1]);
      P->walker[new+j]->potential = e;
      if ( e < P->min_v[l] ){
        P->min_v[l] = e;
        if ( e < P->winner->potential )
          copyBitstring(P->winner, P->walker[new+j]);
      }
      if ( e > P->max_v[l] ) P->max_v[l] = e;
      P->avg_v[l] += e;
      ++(P->psize[l]);
      ++j;
      continue;
    }
    p -= a;
   
    e = mean[P->walker[old+i]->species];
    // Second potential event: walker spawns/dies.
    if ( p < b*(e - P->walker[old+i]->potential) ) {
      copyBitstring(P->walker[new+j], P->walker[old+i]);
      copyBitstring(P->walker[new+j+1], P->walker[old+i]);
      e = P->walker[old+i]->potential;
      if ( e < P->min_v[l] ) P->min_v[l] = e;
      if ( e > P->max_v[l] ) P->max_v[l] = e;
      P->avg_v[l] += 2*e;
      P->psize[l] += 2;
      j += 2;
      continue;
      }
    
    if ( p < b*(P->walker[old+i]->potential - e) ) { // Case of dying.
      continue;
    }
    
    // Third potential event: walker stays.
    copyBitstring(P->walker[new+j], P->walker[old+i]);
    e = P->walker[old+i]->potential;
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
  
  arraysize = 2*nspecies*popsize;
  
  if ( initPopulation(&pop, sat) ) {
    fprintf(stderr,"Could not initialize potential.\n");
    return MEMORY_ERROR;
  }

  seed = time(0);
  //seed = 0 // for testing
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




