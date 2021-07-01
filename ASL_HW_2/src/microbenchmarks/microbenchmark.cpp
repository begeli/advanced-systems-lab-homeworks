#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>

#include "include/microbenchmark.h"
#include "include/tsc_x86.h"
#include "include/foo.h"

#define CYCLES_REQUIRED 1e8
#define REP 40
#define MAX_ELEMENT_VALUE 1000.0
#define ARRAY_SIZE 1000

#define EMPIRICAL_ADDITION_LATENCY 4.0
#define EMPIRICAL_DIVISION_LATENCY 14.0

double* array, *divArr, *fooArr;
double sum, sum0, sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15;
double div0, div1, div2, div3, div4, div5, div6, div7;
double foo0, foo1, foo2, foo3, foo4, foo5, foo6, foo7;
double divRes, fooRes;
bool isCheapDiv;

double generateRandomNo() {
    // random generator from https://stackoverflow.com/questions/55766058/how-can-i-generate-random-doubles-in-c
    srand((unsigned) time(0));
    double randomNo = rand();
    while (randomNo == 0.0) {
      randomNo = rand();
    }
    return (rand() > RAND_MAX / 2 ? -1 : 1) *(MAX_ELEMENT_VALUE / RAND_MAX) * randomNo;
}

void initArray() {
  array = (double*) malloc(ARRAY_SIZE * sizeof(double));
  
  for (int i = 0; i < ARRAY_SIZE; i++) {
    array[i] = generateRandomNo();
  }
}

void initArrayOf1s() {
  array = (double*) malloc(ARRAY_SIZE * sizeof(double));
  
  for (int i = 0; i < ARRAY_SIZE; i++) {
    array[i] = 1.0;
  }
}

void addition_latency() {
    for (int j = 0; j < ARRAY_SIZE; j++) {
      sum += array[j];
    }
}

void divisionLatency() {
  for (int j = 0; j < ARRAY_SIZE; j++) {
    divRes = divRes / array[j];
  }
}

void cheapDivisionLatency() {
  for (int j = 0; j < ARRAY_SIZE; j++) {
      divRes = divRes / array[j];
    }
}

void divisionGap() {
  int j;
  for (j = 0; j < ARRAY_SIZE - 7; j+=8) {
    div0 /= array[j];
    div1 /= array[j+1];
    div2 /= array[j+2];
    div3 /= array[j+3];
    div4 /= array[j+4];
    div5 /= array[j+5];
    div6 /= array[j+6];
    div7 /= array[j+7];
  }
  
  for (; j < ARRAY_SIZE; j++) {
    div0 /= array[j];
  }
}

void minDivisionGap() {
  int j;
  for (j = 0; j < ARRAY_SIZE - 7; j+=8) {
    div0 = div0 / (array[j]);
    div1 = div1 / (array[j+1]);
    div2 = div2 / (array[j+2]);
    div3 = div3 / (array[j+3]);
    div4 = div4 / (array[j+4]);
    div5 = div5 / (array[j+5]);
    div6 = div6 / (array[j+6]);
    div7 = div7 / (array[j+7]);
  }
  
  for (; j < ARRAY_SIZE; j++) {
    div0 = div0 / array[j];
  }
}

void fooLatency() {
  for (int i = 0; i < ARRAY_SIZE; i++) {
    fooRes = foo(fooRes);  
  }
}

void fooGap() {
  int i;
  for (i = 0; i < ARRAY_SIZE - 7; i+=8) {
    foo0 = foo(foo0);
    foo1 = foo(foo1);
    foo2 = foo(foo2);
    foo3 = foo(foo3);
    foo4 = foo(foo4);
    foo5 = foo(foo5);
    foo6 = foo(foo6);
    foo7 = foo(foo7);
  }
  
  for (; i < ARRAY_SIZE; i++) {
    foo0 = foo(foo0);
  }
}


void addition_gap() {
    int j;
    for (j = 0; j < ARRAY_SIZE - 7; j+=8) {
      sum0 += 1;
      sum1 += 1;
      sum2 += 1;
      sum3 += 1;
      sum4 += 1;
      sum5 += 1;
      sum6 += 1;
      sum7 += 1;
    }
    
    for (; j < ARRAY_SIZE; j++) {
      sum0 += 1;
    }
}

void initialize_microbenchmark_data (microbenchmark_mode_t mode) {
  /* You can use to initialize some date if needed */
    switch (mode) {
        case ADD_LAT:
          sum = 0.0;
          initArray();
          break;
        case ADD_GAP:
          //sum = 0.0;
          sum0 = 0.0;
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          sum4 = 0.0;
          sum5 = 0.0;
          sum6 = 0.0;
          sum7 = 0.0;
          initArray();
          break;
        case DIV_LAT:
          divRes = generateRandomNo();
          initArray();
          isCheapDiv = false;
          break;
        case DIV_GAP:
          div0 = generateRandomNo();
          div1 = generateRandomNo();
          div2 = generateRandomNo();
          div3 = generateRandomNo();
          div4 = generateRandomNo();
          div5 = generateRandomNo();
          div6 = generateRandomNo();
          div7 = generateRandomNo();
          isCheapDiv = false;
          initArray();
          break;
        case DIV_LAT_MIN:
          divRes = 0.0;
          initArray();
          isCheapDiv = true;
          break;
        case DIV_GAP_MIN:
          div0 = 0.0;
          div1 = 0.0;
          div2 = 0.0;
          div3 = 0.0;
          div4 = 0.0;
          div5 = 0.0;
          div6 = 0.0;
          div7 = 0.0;
          isCheapDiv = true;
          initArray();
          break;
        case FOO_LAT:
          fooRes = 12321.3123;
          break;
        case FOO_GAP:
          foo0 = generateRandomNo();
          foo1 = generateRandomNo();
          foo2 = generateRandomNo();
          foo3 = generateRandomNo();
          foo4 = generateRandomNo();
          foo5 = generateRandomNo();
          foo6 = generateRandomNo();
          foo7 = generateRandomNo();
          break;
        case FOO_LAT_MIN:
          fooRes = 0.0;
          break;
        case FOO_GAP_MIN:
          foo0 = 0.0;
          foo1 = 0.0;
          foo2 = 0.0;
          foo3 = 0.0;
          foo4 = 0.0;
          foo5 = 0.0;
          foo6 = 0.0;
          foo7 = 0.0;
          break;
        default: break;
    }
}


double microbenchmark_get_add_latency() {
    /* Implement your microbenchmark benchmark here */
    double cycles = 0.0;
    long num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;
    
    // Warm Up
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            addition_latency();
            sum = 0.0;
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
    
    double total_cycles = 0;
    
    for (size_t j = 0; j < REP; j++) {
      start = start_tsc();
      for (int i = 0; i < num_runs; ++i) {
          addition_latency();
          sum = 0.0;
      }
      end = stop_tsc(start);
      cycles = ((double)end) / (num_runs * ARRAY_SIZE);
      total_cycles += cycles;
    }
    
    free(array);
    return total_cycles / REP;
}

double microbenchmark_get_add_gap() {
    /* Implement your microbenchmark benchmark here */
    double cycles = 0.0;
    long num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;
    
    // Warm Up
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            addition_gap();
            sum0 = 0.0;
            sum1 = 0.0;
            sum2 = 0.0;
            sum3 = 0.0;
            /*sum4 = 0.0;
            sum5 = 0.0;
            sum6 = 0.0;
            sum7 = 0.0;*/
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
    
    double total_cycles = 0;
    
    for (size_t j = 0; j < REP; j++) {
      start = start_tsc();
      for (int i = 0; i < num_runs; ++i) {
          addition_gap();
          sum0 = 0.0;
          sum1 = 0.0;
          sum2 = 0.0;
          sum3 = 0.0;
          /* = 0.0;
          sum5 = 0.0;
          sum6 = 0.0;
          sum7 = 0.0;*/
          /*sum8 = 0.0;
          sum9 = 0.0;
          sum10 = 0.0;
          sum11 = 0.0;
          sum12 = 0.0;
          sum13 = 0.0;
          sum14 = 0.0;
          sum15 = 0.0;*/
      }
      end = stop_tsc(start);
      cycles = ((double)end) / (num_runs);
      total_cycles += cycles;
    }
    
    free(array);
    double gap = ((total_cycles / REP) - (EMPIRICAL_ADDITION_LATENCY - 1)) / ARRAY_SIZE;
    return gap; 
  }

double microbenchmark_get_div_latency() {
    /* Implement your microbenchmark benchmark here */
    double cycles = 0.0;
    long num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;
    
    // Warm Up
    if (isCheapDiv) {
      while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            cheapDivisionLatency();
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
      }
    } else {
      while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            divisionLatency();
            divRes = 12345.122;
        }
        cycles = stop_tsc(start);
        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
      }  
    }
    
    
    double total_cycles = 0;
    
    if (isCheapDiv) {
      for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            cheapDivisionLatency();
        }
        end = stop_tsc(start);
        cycles = ((double)end) / (num_runs * ARRAY_SIZE);
        total_cycles += cycles;
      }
    } else {
      for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            divisionLatency();
            divRes = 12345.122;
        }
        end = stop_tsc(start);
        cycles = ((double)end) / (num_runs * ARRAY_SIZE);
        total_cycles += cycles;
      }
    }
    
    free(array);
    return total_cycles / REP;
}

double microbenchmark_get_div_gap() {
    double cycles = 0.0;
    long num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;
    
    // Warm Up
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            divisionGap();
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
    
    double total_cycles = 0;
    
    for (size_t j = 0; j < REP; j++) {
      start = start_tsc();
      for (int i = 0; i < num_runs; ++i) {
          divisionGap();
      }
      end = stop_tsc(start);
      cycles = ((double)end) / (num_runs);
      total_cycles += cycles;
    }
    
    free(array);
    double gap = ((total_cycles / REP) - (EMPIRICAL_DIVISION_LATENCY - 1)) / ARRAY_SIZE;
    return gap;
}

double microbenchmark_get_foo_latency() {
    /* Implement your microbenchmark benchmark here */
    double cycles = 0.0;
    long num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;
    
    // Warm Up
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
          fooLatency();
          //fooRes = 2133113.3213;
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
    
    double total_cycles = 0;
    
    for (size_t j = 0; j < REP; j++) {
      start = start_tsc();
      for (int i = 0; i < num_runs; ++i) {
        fooLatency();
        //fooRes = 2133113.3213;
      }
      end = stop_tsc(start);
      cycles = ((double)end) / (num_runs * ARRAY_SIZE);
      total_cycles += cycles;
    }
    
    //free(array);
    return total_cycles / REP;
}

double microbenchmark_get_foo_gap() {
    /* Implement your microbenchmark benchmark here */
    double cycles = 0.0;
    long num_runs = 1;
    double multiplier = 1;
    myInt64 start, end;
    
    // Warm Up
    while(num_runs < (1 << 14)) {
        start = start_tsc();
        for (int i = 0; i < num_runs; ++i) {
            fooGap();
        }
        cycles = stop_tsc(start);

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }
    
    double total_cycles = 0;
    
    for (size_t j = 0; j < REP; j++) {
      start = start_tsc();
      for (int i = 0; i < num_runs; ++i) {
          fooGap();
      }
      end = stop_tsc(start);
      cycles = ((double)end) / (num_runs);
      total_cycles += cycles;
    }
    
    //free(array);
    double gap = (total_cycles / REP) / ARRAY_SIZE;
    return gap;
}