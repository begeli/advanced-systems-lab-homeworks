#include "common.h"
#include <immintrin.h>
#include <math.h>


void slow_performance1(double *x, double* re, double* im, int N) {
  int k = 0;
  for (int i = 0; i < N; i+=2) {
      re[k] = x[i];
      im[k] = x[i+1];
      k++;
  }
}

void slow_performance2(double *x, double* re, double* im, int N) {
  int k = 0;
  __m256d xv1;
  __m256d xv2;
  __m256d rev;
  __m256d imv;
  
  // Assumed divisibility by 8 which is bad
  for (int i = 0; i < N; i+= 8) {
    xv1 = _mm256_load_pd(x + i);
    xv2 = _mm256_load_pd(x + i + 4);
    
    rev = _mm256_unpacklo_pd(xv1, xv2);
    imv = _mm256_unpackhi_pd(xv1, xv2);
    
    rev = _mm256_permute4x64_pd(rev, 216);
    imv = _mm256_permute4x64_pd(imv, 216);
    
    _mm256_store_pd(re + k, rev);
    _mm256_store_pd(im + k, imv);
    
    k += 4;
  }
}


void maxperformance(double *x, double* re, double* im, int N) {
  int k = 0;
  int i;
  __m256d xv1;
  __m256d xv2;
  __m256d rev;
  __m256d imv;
  
  // kind of unroll the loop 4 times, so... have to do some loop unroll maintanence
  for (i = 0; i < N - 7; i+= 8) { // Also, the code assumes divisibility by 8 which is an incorrect assumption
    xv1 = _mm256_load_pd(x + i);
    xv2 = _mm256_load_pd(x + i + 4);
    
    rev = _mm256_unpacklo_pd(xv1, xv2);
    imv = _mm256_unpackhi_pd(xv1, xv2);
    
    rev = _mm256_permute4x64_pd(rev, 216); // Ugly ass code 11 10 01 00 = 216, also slows down the code a lot, perhaps there is a better way
    imv = _mm256_permute4x64_pd(imv, 216);
    
    _mm256_store_pd(re + k, rev);
    _mm256_store_pd(im + k, imv);
    
    k += 4;
  }
  
  for (; i < N; i += 2) {
    re[k] = x[i];
    im[k] = x[i+1];
    k++;
  }
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() 
{
  add_function(&slow_performance1, "slow_performance1",1);
  add_function(&slow_performance2, "slow_performance2",1);
  add_function(&maxperformance, "maxperformance",1);
}
