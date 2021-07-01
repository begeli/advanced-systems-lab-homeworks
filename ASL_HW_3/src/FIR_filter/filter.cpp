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

void slow_performance2(double *x, double* h, double* y, int N, int M) {
  for (int i = 0; i < N - (M - 1); i++) {
    y[i] = 0.0;
    for (int k = 0; k < M; k++) {
      y[i] += h[k] * fabs(x[i + (M - 1) - k]);
    }
  }
}

void slow_performance3(double *x, double* h, double* y, int N, int M) {
  /* This is the most optimized version. */
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
  //double z[4];
  for (int i = 0; i < N - (M - 1) - 3; i+=4) {
    yVec = _mm256_set1_pd(0.0);
    xVec0 = _mm256_load_pd(x + i);
    xVec0 = _mm256_and_pd(xVec0, maskAbs);
    xVec1 = _mm256_load_pd(x + i + 1);
    xVec1 = _mm256_and_pd(xVec1, maskAbs);
    xVec2 = _mm256_load_pd(x + i + 2);
    xVec2 = _mm256_and_pd(xVec2, maskAbs);
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    
    yVec = _mm256_fmadd_pd( hVec0, xVec3, yVec);
    yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    
    _mm256_store_pd(y + i, yVec);
  }
  //_mm256_store_pd(z, hVec3);
  /*printf("%f %f %f %f\n", z[0], z[1], z[2], z[3]);
  printf("%f %f %f %f\n", h[0], h[1], h[2], h[3]);
  printf("%f %f %f %f %f %f\n", x[0], x[1], x[2], x[3], x[4], x[5]);*/
}

void slow_performance5(double *x, double* h, double* y, int N, int M) {
  __m256d xVec0, xVec1, xVec2, xVec3;
  __m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec;
  
  __m256d xVec4, xVec5, xVec6, xVec7;
  __m256d yVec2;
  
  __m256d xVec8, xVec9, xVec10, xVec11;
  __m256d yVec3;
  
  __m256d xVec12, xVec13, xVec14, xVec15;
  __m256d yVec4;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  //double z[4];
  int i;
  for ( i = 0; i < N - (M - 1) - 15; i+=16) {
    //yVec = _mm256_set1_pd(0.0);
    //yVec2 = _mm256_set1_pd(0.0);
    //yVec3 = _mm256_set1_pd(0.0);
    //yVec4 = _mm256_set1_pd(0.0);
    
    xVec0 = _mm256_load_pd(x + i);
    xVec0 = _mm256_and_pd(xVec0, maskAbs);
    xVec1 = _mm256_load_pd(x + i + 1);
    xVec1 = _mm256_and_pd(xVec1, maskAbs);
    xVec2 = _mm256_load_pd(x + i + 2);
    xVec2 = _mm256_and_pd(xVec2, maskAbs);
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    
    xVec4 = _mm256_load_pd(x + i + 4);
    xVec4 = _mm256_and_pd(xVec4, maskAbs);
    xVec5 = _mm256_load_pd(x + i + 5);
    xVec5 = _mm256_and_pd(xVec5, maskAbs);
    xVec6 = _mm256_load_pd(x + i + 6);
    xVec6 = _mm256_and_pd(xVec6, maskAbs);
    xVec7 = _mm256_load_pd(x + i + 7);
    xVec7 = _mm256_and_pd(xVec7, maskAbs);
    
    xVec8 = _mm256_load_pd(x + i + 8);
    xVec8 = _mm256_and_pd(xVec8, maskAbs);
    xVec9 = _mm256_load_pd(x + i + 9);
    xVec9 = _mm256_and_pd(xVec9, maskAbs);
    xVec10 = _mm256_load_pd(x + i + 10);
    xVec10 = _mm256_and_pd(xVec10, maskAbs);
    xVec11 = _mm256_load_pd(x + i + 11);
    xVec11 = _mm256_and_pd(xVec11, maskAbs);
    
    xVec12 = _mm256_load_pd(x + i + 12);
    xVec12 = _mm256_and_pd(xVec12, maskAbs);
    xVec13 = _mm256_load_pd(x + i + 13);
    xVec13 = _mm256_and_pd(xVec13, maskAbs);
    xVec14 = _mm256_load_pd(x + i + 14);
    xVec14 = _mm256_and_pd(xVec14, maskAbs);
    xVec15 = _mm256_load_pd(x + i + 15);
    xVec15 = _mm256_and_pd(xVec15, maskAbs);
    
    yVec = _mm256_mul_pd( hVec0, xVec3);
    yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    
    yVec2 = _mm256_mul_pd( hVec0, xVec7);
    yVec2 = _mm256_fmadd_pd( hVec1, xVec6, yVec2);
    yVec2 = _mm256_fmadd_pd( hVec2, xVec5, yVec2);
    yVec2 = _mm256_fmadd_pd( hVec3, xVec4, yVec2);
    
    yVec3 = _mm256_mul_pd( hVec0, xVec11);
    yVec3 = _mm256_fmadd_pd( hVec1, xVec10, yVec3);
    yVec3 = _mm256_fmadd_pd( hVec2, xVec9, yVec3);
    yVec3 = _mm256_fmadd_pd( hVec3, xVec8, yVec3);
    
    yVec4 = _mm256_mul_pd( hVec0, xVec15);
    yVec4 = _mm256_fmadd_pd( hVec1, xVec14, yVec4);
    yVec4 = _mm256_fmadd_pd( hVec2, xVec13, yVec4);
    yVec4 = _mm256_fmadd_pd( hVec3, xVec12, yVec4);
    
    _mm256_store_pd(y + i, yVec);
    _mm256_store_pd(y + i + 4, yVec2);
    _mm256_store_pd(y + i + 8, yVec3);
    _mm256_store_pd(y + i + 12, yVec4);
  }
  
  for (; i < N - (M - 1); i++) {
    y[i] = 0.0;
    for (int k = 0; k < M; k++) {
      y[i] += h[k] * fabs(x[i + (M - 1) - k]);
    }
  }
}

void slow_performance6(double *x, double* h, double* y, int N, int M) {
  __m256d xVec0, xVec1, xVec2, xVec3;
  __m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec;
  
  __m256d xVec4, xVec5, xVec6, xVec7;
  //__m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec2;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  __m256d zeros = _mm256_setzero_pd();
  //double z[4];
  int i;
  for (i = 0; i < N - (M - 1) - 7; i+=8) {
    //yVec = _mm256_setzero_pd();//_mm256_set1_pd(0.0);
    //yVec2 = _mm256_setzero_pd();//_mm256_set1_pd(0.0);
    
    xVec0 = _mm256_load_pd(x + i);
    xVec0 = _mm256_and_pd(xVec0, maskAbs);
    xVec4 = _mm256_load_pd(x + i + 4);
    xVec4 = _mm256_and_pd(xVec4, maskAbs);
    xVec7 = _mm256_load_pd(x + i + 7);
    xVec7 = _mm256_and_pd(xVec7, maskAbs);
    
    /*xVec2 = _mm256_load_pd(x + i + 2);
    xVec2 = _mm256_and_pd(xVec2, maskAbs);*/
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    xVec5 = _mm256_load_pd(x + i + 5);
    xVec5 = _mm256_and_pd(xVec5, maskAbs);
    xVec6 = _mm256_load_pd(x + i + 6);
    xVec6 = _mm256_and_pd(xVec6, maskAbs);
    
    
    xVec1 = _mm256_blend_pd(xVec0, xVec4, 1);
    xVec1 = _mm256_permute4x64_pd(xVec1, 57);
    xVec2 = _mm256_blend_pd(xVec0, xVec4, 3);
    xVec2 = _mm256_permute4x64_pd(xVec2, 78);
    
    //yVec = _mm256_fmadd_pd( hVec0, xVec3, yVec);
    yVec = _mm256_fmadd_pd( hVec0, xVec3, zeros);
    yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    
    //yVec2 = _mm256_fmadd_pd( hVec0, xVec7, yVec2);
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
  //__m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec2;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  //__m256d zeros = _mm256_setzero_pd();
  //double z[4];
  int i;
  for (i = 0; i < N - (M - 1) - 7; i+=8) {
    //yVec = _mm256_setzero_pd();//_mm256_set1_pd(0.0);
    //yVec2 = _mm256_setzero_pd();//_mm256_set1_pd(0.0);
    
    xVec0 = _mm256_load_pd(x + i);
    //xVec0 = _mm256_and_pd( xVec0, maskAbs);
    yVec = _mm256_mul_pd( hVec3, _mm256_and_pd( xVec0, maskAbs));
    
    xVec1 = _mm256_load_pd(x + i + 1);
    //xVec1 = _mm256_and_pd(xVec1, maskAbs);
    yVec = _mm256_fmadd_pd( hVec2, _mm256_and_pd(xVec1, maskAbs), yVec);
    
    xVec2 = _mm256_load_pd(x + i + 2);
    //xVec2 = _mm256_and_pd(xVec2, maskAbs);
    yVec = _mm256_fmadd_pd( hVec1, _mm256_and_pd(xVec2, maskAbs), yVec);
    
    xVec3 = _mm256_load_pd(x + i + 3);
    xVec3 = _mm256_and_pd(xVec3, maskAbs);
    yVec = _mm256_fmadd_pd( hVec0, xVec3, yVec);
    
    xVec4 = _mm256_load_pd(x + i + 4);
    xVec4 = _mm256_and_pd(xVec4, maskAbs);
    xVec5 = _mm256_load_pd(x + i + 5);
    //xVec5 = _mm256_and_pd(xVec5, maskAbs);
    xVec6 = _mm256_load_pd(x + i + 6);
    xVec6 = _mm256_and_pd(xVec6, maskAbs);
    xVec7 = _mm256_load_pd(x + i + 7);
    xVec7 = _mm256_and_pd(xVec7, maskAbs);
    
    //yVec = _mm256_fmadd_pd( hVec0, xVec3, yVec);
    //yVec = _mm256_mul_pd( hVec0, xVec3);
    //yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    //yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    //yVec = _mm256_fmadd_pd( hVec2, _mm256_and_pd(xVec1, maskAbs), yVec);
    //yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    //yVec = _mm256_fmadd_pd( hVec3, _mm256_and_pd(xVec0, maskAbs), yVec);
    
    //yVec2 = _mm256_fmadd_pd( hVec0, xVec7, yVec2);
    yVec2 = _mm256_mul_pd( hVec0, xVec7);
    yVec2 = _mm256_fmadd_pd( hVec1, xVec6, yVec2);
    //yVec2 = _mm256_fmadd_pd( hVec2, xVec5, yVec2);
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
  //__m256d hVec0, hVec1, hVec2, hVec3;
  __m256d yVec2;
  
  hVec0 = _mm256_set1_pd(h[0]);
  hVec1 = _mm256_set1_pd(h[1]);
  hVec2 = _mm256_set1_pd(h[2]);
  hVec3 = _mm256_set1_pd(h[3]);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  //__m256d zeros = _mm256_setzero_pd();
  //double z[4];
  int i;
  for (i = 0; i < N - (M - 1) - 7; i+=8) {
    //yVec = _mm256_setzero_pd();//_mm256_set1_pd(0.0);
    //yVec2 = _mm256_setzero_pd();//_mm256_set1_pd(0.0);
    
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
    //xVec5 = _mm256_and_pd(xVec5, maskAbs);
    xVec6 = _mm256_load_pd(x + i + 6);
    xVec6 = _mm256_and_pd(xVec6, maskAbs);
    xVec7 = _mm256_load_pd(x + i + 7);
    xVec7 = _mm256_and_pd(xVec7, maskAbs);
    
    //yVec = _mm256_fmadd_pd( hVec0, xVec3, yVec);
    yVec = _mm256_mul_pd( hVec0, xVec3);
    yVec = _mm256_fmadd_pd( hVec1, xVec2, yVec);
    yVec = _mm256_fmadd_pd( hVec2, xVec1, yVec);
    //yVec = _mm256_fmadd_pd( hVec2, _mm256_and_pd(xVec1, maskAbs), yVec);
    yVec = _mm256_fmadd_pd( hVec3, xVec0, yVec);
    //yVec = _mm256_fmadd_pd( hVec3, _mm256_and_pd(xVec0, maskAbs), yVec);
    
    //yVec2 = _mm256_fmadd_pd( hVec0, xVec7, yVec2);
    yVec2 = _mm256_mul_pd( hVec0, xVec7);
    yVec2 = _mm256_fmadd_pd( hVec1, xVec6, yVec2);
    //yVec2 = _mm256_fmadd_pd( hVec2, xVec5, yVec2);
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
  add_function(&slow_performance2, "slow_performance2",1);
  add_function(&slow_performance4, "slow_performance4",1); // vectorized
  //add_function(&slow_performance5, "slow_performance5",1); //unrollimg more than 2 times slows it down
  //add_function(&slow_performance3, "slow_performance3",1);
  add_function(&slow_performance6, "slow_performance6",1);
  add_function(&slow_performance7, "slow_performance7",1);
  add_function(&maxperformance, "maxperformance",1); // unrolling 2 times
}