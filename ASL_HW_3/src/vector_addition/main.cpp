/**
*      _________   _____________________  ____  ______
*     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
*    / /_  / /| | \__ \ / / / /   / / / / / / / __/
*   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
*  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
*
*  http://www.acl.inf.ethz.ch/teaching/fastcode
*  How to Write Fast Numerical Code 263-2300 - ETH Zurich
*  Copyright (C) 2019 
*                   Tyler Smith        (smitht@inf.ethz.ch) 
*                   Alen Stojanov      (astojanov@inf.ethz.ch)
*                   Gagandeep Singh    (gsingh@inf.ethz.ch)
*                   Markus Pueschel    (pueschel@inf.ethz.ch)
*
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program. If not, see http://www.gnu.org/licenses/.
*/

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <random>
#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "tsc_x86.h"

using namespace std;

#define NR 32
#define CYCLES_REQUIRED (1*1e8)
#define REP 50
#define EPS (1e-3)

void kernel_base(vector_t x, vector_t y, vector_t z, int n) {
  for (int i = 0; i < n; i++) {
    if (x.id[i] == y.id[i]) {
      z.val[i] = x.val[i] + y.val[i];
      z.id[i]  = x.id[i];
    }
    else if (fabs(x.val[i]) > fabs(y.val[i])){
      z.val[i]  = x.val[i];
      z.id[i]   = x.id[i]; 
      z.val[n] += fabs(y.val[i]);
    }
    else {
      z.val[i]  = y.val[i];
      z.id[i]   = y.id[i]; 
      z.val[n] += fabs(x.val[i]);
    }
  }
}

/* prototype of the function you need to optimize */
// typedef void(*comp_func)(vector_t, vector_t, vector_t, int);

//headers
void   register_functions();
double get_perf_score(comp_func f);
double perf_test(comp_func f, string desc, int flops);
void add_function(comp_func f, string name, int flop);

/* Global vars, used to keep track of student functions */
vector<comp_func> userFuncs;
vector<string> funcNames;
vector<int> funcFlops;
int numFuncs = 0;


template<typename T>
void rands(T * m, size_t row, size_t col)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<T> dist(-1.0, 1.0);
    for (size_t i = 0; i < row*col; ++i)  
        m[i] = dist(gen);
}

void rands_id(uint64_t * m, size_t n)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_int_distribution<uint64_t> dist(0, n-1);
    for (size_t i = 0; i < n; ++i) {
      m[i] = i;
    }
    
    for (size_t i = 0; i < 10*n / 16; ++i) {
      size_t j = dist(gen);
      m[j]     = dist(gen);
    }
}

void rands_vec(vector_t& v, size_t n)
{
    rands(v.val, 1, n);
    rands_id(v.id, n);
}

template<typename T>
void build(T **a, int m, int n)
{
    *a = static_cast<T *>(aligned_alloc(32, m * n * sizeof(T)));
    rands(*a, m, n);
}

void build_vec(vector_t& v, int n)
{
    v.val = static_cast<double *>  (aligned_alloc(32, n * sizeof(double)));
    v.id  = static_cast<uint64_t *>(aligned_alloc(32, n * sizeof(uint64_t)));
    rands_vec(v, n);
}

template<typename T>
void destroy(T* m)
{
    free(m);
}

void destroy_vec(vector_t& m)
{
    free(m.val);
    free(m.id);
}

template<typename T>
T nrm_sqr_diff(T *x, T *y, int n) {
    T nrm_sqr = 0.0;
    for(int i = 0; i < n; i++) {
        nrm_sqr += (x[i] - y[i]) * (x[i] - y[i]);
    }
    
    if (isnan(nrm_sqr)) {
      nrm_sqr = INFINITY;
    }
    
    return nrm_sqr;
}


double nrm_sqr_diff_vec(vector_t& x, vector_t& y, int n) {
    double nrm_sqr = 0.0;
    for(int i = 0; i < n; i++) {
        nrm_sqr += (x.val[i] - y.val[i]) * (x.val[i] - y.val[i]);
        
        if (x.id[i] != y.id[i]) {
          nrm_sqr = INFINITY;
        }
    }
    
    if (isnan(nrm_sqr)) {
      nrm_sqr = INFINITY;
    }
    
    return nrm_sqr;
}

/*
* Registers a user function to be tested by the driver program. Registers a
* string description of the function as well
*/
void add_function(comp_func f, string name, int flops)
{
    userFuncs.push_back(f);
    funcNames.emplace_back(name);
    funcFlops.push_back(flops);

    numFuncs++;
}


/*
* Checks the given function for validity. If valid, then computes and
* reports and returns the number of cycles required per iteration
*/
double perf_test(comp_func f, string desc, int flops)
{
    double cycles = 0.;
    long num_runs = 100;
    double multiplier = 1;
    myInt64 start, end;
    
    const int n = 512;
    static double   x_val[n]   __attribute__ ((aligned (32)));
    static double   y_val[n]   __attribute__ ((aligned (32)));
    static double   z_val[n+1] __attribute__ ((aligned (32)));
    static uint64_t x_id[n]    __attribute__ ((aligned (32)));
    static uint64_t y_id[n]    __attribute__ ((aligned (32)));
    static uint64_t z_id[n+1]  __attribute__ ((aligned (32)));
    vector_t x = { .val = x_val, .id = x_id };
    vector_t y = { .val = y_val, .id = y_id };
    vector_t z = { .val = z_val, .id = z_id };
    rands_vec(x, n);
    rands_vec(y, n);
    rands_vec(z, n+1);

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            f(x, y, z, n); 
            x.val[0] = z.val[n];          
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);
        
    } while (multiplier > 2);
  

    list<double> cyclesList;

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {

        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            f(x, y, z, n);
            x.val[0] = z.val[n];
        }
        end = stop_tsc(start);
        
        /* Try other values */
        rands_vec(x, n);
        rands_vec(y, n);

        cycles = ((double)end) / num_runs;
        total_cycles += cycles;

        cyclesList.push_back(cycles);
    }
    total_cycles /= REP;

    cycles = total_cycles;
    cyclesList.sort();
    return  cycles; //cyclesList.front();
}

int main(int argc, char **argv)
{
  cout << "Starting program. ";
  double perf;
  int i;

  register_functions();

  if (numFuncs == 0){
    cout << endl;
    cout << "No functions registered - nothing for driver to do" << endl;
    cout << "Register functions by calling register_func(f, name)" << endl;
    cout << "in register_funcs()" << endl;

    return 0;
  }
  cout << numFuncs << " functions registered." << endl;
   
    //Check validity of functions. 
  int n_base = 1024;
  vector_t x, y, z, z_base, z_old;
  build_vec(x, n_base);
  build_vec(y, n_base);
  build_vec(z, n_base+1);
  build_vec(z_base, n_base+1);
  build_vec(z_old,  n_base+1);
  
  for (int k = 0; k < 1; k++) {
    int n = n_base - k;
    memcpy(z_old.val,  z.val, (n+1)*sizeof(double));
    memcpy(z_old.id,   z.id,  (n+1)*sizeof(uint64_t));
    kernel_base(x, y, z, n);
    memcpy(z_base.val, z.val, (n+1)*sizeof(double));
    memcpy(z_base.id,  z.id,  (n+1)*sizeof(uint64_t));
  
    for (i = 0; i < numFuncs; i++) {
      memcpy(z.val, z_old.val, (n+1)*sizeof(double));
      memcpy(z.id , z_old.id,  (n+1)*sizeof(uint64_t));
      comp_func f = userFuncs[i];
      f(x, y, z, n);
      double error = nrm_sqr_diff_vec(z, z_base, n+1);
      
      if (error > EPS) {
        cout << error << endl;
        cout << "ERROR!!!!  the results for the " << i+1 << "th function are different to the previous" << std::endl;
        // return 0;
      }
    }
  }
  destroy_vec(x);
  destroy_vec(y);
  destroy_vec(z);
  destroy_vec(z_base);
  destroy_vec(z_old);


  for (i = 0; i < numFuncs; i++) {
    perf = perf_test(userFuncs[i], funcNames[i], 1);
    cout << endl << "Running: " << funcNames[i] << endl;
    cout << perf << " cycles" << endl;
  }

  return 0;
}