#include "common.h"
#include <immintrin.h>
#include <math.h>

void slow_performance1(double *x, double* h, double* y, int N, int M) {
  for (int i = 0; i < N - (M - 1); i++) {
    y[i] = 0.0;
    for (int k = 0; k < M; k++) {
      y[i] += h[k] * fabs(x[i + (M - 1) - k]);
    }
  }
}

void slow_performance4(double *x, double* h, double* y, int N, int M) {
  __m256d xVec0, xVec1, xVec2, xVec3;
  __m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  for (int i = 0; i < N - (M - 1) - 3; i+=4) {
    yVec = _mm256_set1_pd(0.0);
    
    // Load vectors and take their absolute value
    xVec0 = _mm256_load_pd(x + i);
    xVec0 = _mm256_and_pd(xVec0, maskAbs);
    xVec1 = _mm256_load_pd(x + i + 1);
    xVec1 = _mm256_and_pd(xVec1, maskAbs);
    xVec2 = _mm256_load_pd(x + i + 2);
    xVec2 = _mm256_and_pd(xVec2, maskAbs);
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    
    // Compute the final value of the output
    yVec = _mm256_fmadd_pd( hVec0, xVec3, yVec);
    yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    
    _mm256_store_pd(y + i, yVec);
  }
}

void slow_performance6(double *x, double* h, double* y, int N, int M) {
  __m256d xVec0, xVec1, xVec2, xVec3;
  __m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec;
  
  __m256d xVec4, xVec5, xVec6, xVec7;
  __m256d yVec2;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  __m256d zeros = _mm256_setzero_pd();
  
  int i;
  for (i = 0; i < N - (M - 1) - 7; i+=8) {
    // Load vectors and take their absolute value
    xVec0 = _mm256_load_pd(x + i);
    xVec0 = _mm256_and_pd(xVec0, maskAbs);
    xVec4 = _mm256_load_pd(x + i + 4);
    xVec4 = _mm256_and_pd(xVec4, maskAbs);
    xVec7 = _mm256_load_pd(x + i + 7);
    xVec7 = _mm256_and_pd(xVec7, maskAbs);
    
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    xVec5 = _mm256_load_pd(x + i + 5);
    xVec5 = _mm256_and_pd(xVec5, maskAbs);
    xVec6 = _mm256_load_pd(x + i + 6);
    xVec6 = _mm256_and_pd(xVec6, maskAbs);
    
    // No idea what I am doin here... But judging by the fact that I did not do this in the other versions suggest that this didn't work.
    xVec1 = _mm256_blend_pd(xVec0, xVec4, 1);
    xVec1 = _mm256_permute4x64_pd(xVec1, 57);
    xVec2 = _mm256_blend_pd(xVec0, xVec4, 3);
    xVec2 = _mm256_permute4x64_pd(xVec2, 78);
    
    // Compute the final value of the output
    yVec = _mm256_fmadd_pd( hVec0, xVec3, zeros);
    yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    
    yVec2 = _mm256_fmadd_pd( hVec0, xVec7, zeros);
    yVec2 = _mm256_fmadd_pd( hVec1, xVec6, yVec2);
    yVec2 = _mm256_fmadd_pd( hVec2, xVec5, yVec2);
    yVec2 = _mm256_fmadd_pd( hVec3, xVec4, yVec2);
    
    _mm256_store_pd(y + i, yVec);
    _mm256_store_pd(y + i + 4, yVec2);
  }
  
  for (; i < N - (M - 1); i++) {
    y[i] = 0.0;
    for (int k = 0; k < M; k++) {
      y[i] += h[k] * fabs(x[i + (M - 1) - k]);
    }
  }
}

void slow_performance7(double *x, double* h, double* y, int N, int M) {
  __m256d xVec0, xVec1, xVec2, xVec3;
  __m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec;
  
  __m256d xVec4, xVec5, xVec6, xVec7;
  __m256d yVec2;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  int i;
  for (i = 0; i < N - (M - 1) - 7; i+=8) {
    // Load vectors and take their absolute value / Also start calculating the the output the interleave instructions
    xVec0 = _mm256_load_pd(x + i);
    yVec = _mm256_mul_pd( hVec3, _mm256_and_pd( xVec0, maskAbs));
    
    xVec1 = _mm256_load_pd(x + i + 1);
    yVec = _mm256_fmadd_pd( hVec2, _mm256_and_pd(xVec1, maskAbs), yVec);
    
    xVec2 = _mm256_load_pd(x + i + 2);
    yVec = _mm256_fmadd_pd( hVec1, _mm256_and_pd(xVec2, maskAbs), yVec);
    
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    yVec = _mm256_fmadd_pd( hVec0, xVec3, yVec);
    
    // Load vectors and take their absolute value 
    xVec4 = _mm256_load_pd(x + i + 4);
    xVec4 = _mm256_and_pd(xVec4, maskAbs);
    xVec5 = _mm256_load_pd(x + i + 5);
    xVec6 = _mm256_load_pd(x + i + 6);
    xVec6 = _mm256_and_pd(xVec6, maskAbs);
    xVec7 = _mm256_load_pd(x + i + 7);
    xVec7 = _mm256_and_pd(xVec7, maskAbs);
    
    // Compute the final value of the output
    yVec2 = _mm256_mul_pd( hVec0, xVec7);
    yVec2 = _mm256_fmadd_pd( hVec1, xVec6, yVec2);
    yVec2 = _mm256_fmadd_pd( hVec2, _mm256_and_pd(xVec5, maskAbs), yVec2);
    yVec2 = _mm256_fmadd_pd( hVec3, xVec4, yVec2);
    
    _mm256_store_pd(y + i, yVec);
    _mm256_store_pd(y + i + 4, yVec2);
  }
  
  for (; i < N - (M - 1); i++) {
    y[i] = 0.0;
    for (int k = 0; k < M; k++) {
      y[i] += h[k] * fabs(x[i + (M - 1) - k]);
    }
  }
}

void maxperformance(double *x, double* h, double* y, int N, int M) {
  __m256d xVec0, xVec1, xVec2, xVec3;
  __m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec;
  
  __m256d xVec4, xVec5, xVec6, xVec7;
  __m256d yVec2;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  int i;
  for (i = 0; i < N - (M - 1) - 7; i+=8) {
    // Load vectors and take their absolute value 
    xVec0 = _mm256_load_pd(x + i);
    xVec0 = _mm256_and_pd( xVec0, maskAbs);
    xVec1 = _mm256_load_pd(x + i + 1);
    xVec1 = _mm256_and_pd(xVec1, maskAbs);
    xVec2 = _mm256_load_pd(x + i + 2);
    xVec2 = _mm256_and_pd(xVec2, maskAbs);
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    
    xVec4 = _mm256_load_pd(x + i + 4);
    xVec4 = _mm256_and_pd(xVec4, maskAbs);
    xVec5 = _mm256_load_pd(x + i + 5);
    xVec6 = _mm256_load_pd(x + i + 6);
    xVec6 = _mm256_and_pd(xVec6, maskAbs);
    xVec7 = _mm256_load_pd(x + i + 7);
    xVec7 = _mm256_and_pd(xVec7, maskAbs);
    
    // Compute the final value of the output
    yVec = _mm256_mul_pd( hVec0, xVec3);
    yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    
    yVec2 = _mm256_mul_pd( hVec0, xVec7);
    yVec2 = _mm256_fmadd_pd( hVec1, xVec6, yVec2);
    yVec2 = _mm256_fmadd_pd( hVec2, _mm256_and_pd(xVec5, maskAbs), yVec2);
    yVec2 = _mm256_fmadd_pd( hVec3, xVec4, yVec2);
    
    _mm256_store_pd(y + i, yVec);
    _mm256_store_pd(y + i + 4, yVec2);
  }
  
  for (; i < N - (M - 1); i++) {
    y[i] = 0.0;
    for (int k = 0; k < M; k++) {
      y[i] += h[k] * fabs(x[i + (M - 1) - k]);
    }
  }
}

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() 
{
  add_function(&slow_performance1, "slow_performance1",1);
  add_function(&slow_performance4, "slow_performance4",1); // vectorized
  add_function(&slow_performance6, "slow_performance6",1);
  add_function(&slow_performance7, "slow_performance7",1);
  add_function(&maxperformance, "maxperformance",1); // unrolling 2 times
}
