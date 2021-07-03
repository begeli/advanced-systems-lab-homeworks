#include "common.h"
#include <math.h>
#include <iostream>

#define C1 0.2
#define C2 0.3

void slow_performance1(float* x, float *y, float *z, int n) {
    for (int i = 0; i < n - 2; i++) {
        x[i]     = x[i] / M_SQRT2 + y[i] * C1;
        x[i+1]  += z[(i % 4) * 10] * C2;
        x[i+2]  += sin ((2 * M_PI * i) / 3) * y[i+2];
    }
}

void slow_performance2(float* x, float *y, float *z, int n) {
  /* Scalar replacement */
  float yi, yip2, xi, xip1, xip2, zi; 
  for (int i = 0; i < n - 2; i++) {
        // Load
        yi = y[i];
        yip2 = y[i+2];
        xi = x[i];
        xip1 = x[i+1];
        xip2 = x[i+2];
        zi = z[(i % 4) * 10];
        
        // Compute
        xi = xi / M_SQRT2 + yi * C1;
        xip1 = xip1 + zi * C2;
        xip2 = xip2 + sin ((2 * M_PI * i) / 3) * yip2;
        
        // Store
        x[i]  = xi;
        x[i+1]  = xip1;
        x[i+2]  = xip2;
    }
}

void slow_performance3(float* x, float *y, float *z, int n) {
  //Scalar replacement & Loop Unrolling by a factor of hmmm lets say 6
  float yi, yip2, xi, xip1, xip2, zi; 
  float xip3, xip4, xip5, xip6, xip7;
  float yip1, yip3, yip4, yip5, yip6, yip7;
  float zi1, zi2, zi3;
  
  float t;
  
  const int UNROLL_AMOUNT = 6;
  int limit = n - (UNROLL_AMOUNT - 1);
  int i;
  for (i = 0; i < limit - 2; i+=UNROLL_AMOUNT) {
      // Load 
      yi = y[i];
      yip1 = y[i+1];
      yip2 = y[i+2];
      yip3 = y[i+3];
      yip4 = y[i+4];
      yip5 = y[i+5];
      yip6 = y[i+6];
      yip7 = y[i+7];
      zi = z[(i % 4) * 10];
      zi1 = z[((i + 1) % 4) * 10];
      zi2 = z[((i + 2) % 4) * 10];
      zi3 = z[((i + 3) % 4) * 10];
      xi = x[i];
      xip1 = x[i+1];
      xip2 = x[i+2];
      xip3 = x[i+3];
      xip4 = x[i+4];
      xip5 = x[i+5];
      xip6 = x[i+6];
      xip7 = x[i+7];
      
      // Compute
      t = yi * C1;
      xi = xi / M_SQRT2;
      xi = xi + t;
      t = zi * C2;
      xip1 = xip1 + t;
      t = sin ((2 * M_PI * i) / 3) * yip2;
      xip2  = xip2 + t;
      
      t = yip1 * C1;
      xip1 = xip1 / M_SQRT2;
      xip1 = xip1 + t;
      t = zi1 * C2;
      xip2  = xip2 + t;
      t = sin ((2 * M_PI * (i+1)) / 3) * yip3;
      xip3 = xip3 + t;
      
      t = yip2 * C1;
      xip2 = xip2 / M_SQRT2;
      xip2 = xip2 + t;
      t = zi2 * C2;
      xip3  = xip3 + t;
      t = sin ((2 * M_PI * (i+2)) / 3) * yip4;
      xip4 = xip4 + t;
      
      t = yip3 * C1;
      xip3 = xip3 / M_SQRT2;
      xip3  = xip3 + t;
      t = zi3 * C2;
      xip4  = xip4 + t;
      t = sin ((2 * M_PI * (i+3)) / 3) * yip5;
      xip5  = xip5 + sin ((2 * M_PI * (i+3)) / 3) * yip5;
      
      t = yip4 * C1;
      xip4 = xip4 / M_SQRT2;
      xip4 = xip4 + t;
      t = zi * C2;
      xip5 = xip5 + t;
      t = sin ((2 * M_PI * (i+4)) / 3) * yip6;
      xip6 = xip6 + t;
      
      t = yip5 * C1;
      xip5 = xip5 / M_SQRT2;
      xip5 = xip5 + t;
      t = zi1 * C2;
      xip6 = xip6 + t;
      t = sin ((2 * M_PI * (i+5)) / 3) * yip7;
      xip7 = xip7 + t;
      
      // Store
      x[i] = xi;
      x[i+1] = xip1;
      x[i+2] = xip2;
      x[i+3] = xip3;
      x[i+4] = xip4;
      x[i+5] = xip5;
      x[i+6] = xip6;
      x[i+7] = xip7;
  }
  
  for (; i < n - 2; i++) {
      // Load
      yi = y[i];
      yip2 = y[i+2];
      xi = x[i];
      xip1 = x[i+1];
      xip2 = x[i+2];
      zi = z[(i % 4) * 10]; // take outside of loop and calculate z[0], z[10], z[20], z[30]
      
      // Compute
      xi = xi / M_SQRT2 + yi * C1;
      xip1 = xip1 + zi * C2;
      xip2 = xip2 + sin ((2 * M_PI * i) / 3) * yip2; // take sin out and calculate the possible values 0, 60, 120 etc...
      
      // Store
      x[i]  = xi;//x[i] / M_SQRT2 + yi * C1;
      x[i+1]  = xip1;
      x[i+2]  = xip2;
  }
}

void slow_performance4(float* x, float *y, float *z, int n) {
  /* Factor out common expression */
  double mt_mpi_dt = 2 * M_PI / 3;
  for (int i = 0; i < n - 2; i++) {
      x[i]  = x[i] / M_SQRT2 + y[i] * C1; // removed divison here
      x[i+1]  += z[(i % 4) * 10] * C2;
      x[i+2]  += sin ( mt_mpi_dt * i) * y[i+2];
  }
}


void slow_performance6(float* x, float *y, float *z, int n) {
  //Scalar replacement & Loop Unrolling by a factor of hmmm lets say 6
  float yi, yip2, xi, xip1, xip2; 
  float xip3, xip4, xip5, xip6, xip7;
  float yip1, yip3, yip4, yip5, yip6, yip7;
  float zi, zi1, zi2, zi3;
  double sin0, sin1, sin2, sin3, sin4, sin5, sin6;
  zi = z[0];
  zi1 = z[10];
  zi2 = z[20];
  zi3 = z[30];
  float t;
  
  double mt_mpi_dt = 2 * M_PI / 3;
  double* sineValues = (double*) malloc(6 * sizeof(double));
  sineValues[0] = sin (0);
  sineValues[1] = sin (mt_mpi_dt);
  sineValues[2] = sin (mt_mpi_dt * 2);
  sineValues[3] = sin (mt_mpi_dt * 3);
  sineValues[4] = sin (mt_mpi_dt * 4);
  sineValues[5] = sin (mt_mpi_dt * 5);
  
  const int UNROLL_AMOUNT = 6;
  int limit = n - (UNROLL_AMOUNT - 1);
  int i;
  for (i = 0; i < limit - 2; i+=UNROLL_AMOUNT) {
      // Load 
      yi = y[i];
      yip1 = y[i+1];
      yip2 = y[i+2];
      yip3 = y[i+3];
      yip4 = y[i+4];
      yip5 = y[i+5];
      yip6 = y[i+6];
      yip7 = y[i+7];
      zi = z[(i % 4) * 10];
      zi1 = z[((i + 1) % 4) * 10];
      zi2 = z[((i + 2) % 4) * 10];
      zi3 = z[((i + 3) % 4) * 10];
      xi = x[i];
      xip1 = x[i+1];
      xip2 = x[i+2];
      xip3 = x[i+3];
      xip4 = x[i+4];
      xip5 = x[i+5];
      xip6 = x[i+6];
      xip7 = x[i+7];
      sin0 = sineValues[0];
      sin1 = sineValues[1];
      sin2 = sineValues[2];
      sin3 = sineValues[3];
      sin4 = sineValues[4];
      sin5 = sineValues[5];
      
      // Compute
      t = yi * 0.2;//C1;
      xi = xi * M_SQRT1_2;
      xi = xi + t;
      t = zi * 0.3;
      xip1 = xip1 + t;
      t = sin0 * yip2;//sin (mt_mpi_dt * i ) * yip2;
      xip2  = xip2 + t;
      
      t = yip1 * 0.2;
      xip1 = xip1 * M_SQRT1_2;
      xip1 = xip1 + t;
      t = zi1 * 0.3;
      xip2  = xip2 + t;
      t = sin1 * yip3;
      xip3 = xip3 + t;
      
      t = yip2 * 0.2;
      xip2 = xip2 * M_SQRT1_2;
      xip2 = xip2 + t;
      t = zi2 * 0.3;
      xip3  = xip3 + t;
      t = sin2 * yip4;
      xip4 = xip4 + t;
      
      t = yip3 * 0.2;
      xip3 = xip3 * M_SQRT1_2;
      xip3  = xip3 + t;
      t = zi3 * 0.3;
      xip4  = xip4 + t;
      t = sin3 * yip5;
      xip5  = xip5 + t;
      
      t = yip4 * 0.2;
      xip4 = xip4 * M_SQRT1_2;
      xip4 = xip4 + t;
      t = zi * 0.3;
      xip5 = xip5 + t;
      t = sin4 * yip6;
      xip6 = xip6 + t;
      
      t = yip5 * 0.2;
      xip5 = xip5 * M_SQRT1_2;
      xip5 = xip5 + t;
      t = zi1 * 0.3;
      xip6 = xip6 + t;
      t = sin5 * yip7;
      xip7 = xip7 + t;
      
      // Store
      x[i] = xi;
      x[i+1] = xip1;
      x[i+2] = xip2;
      x[i+3] = xip3;
      x[i+4] = xip4;
      x[i+5] = xip5;
      x[i+6] = xip6;
      x[i+7] = xip7;
  }
  
  for (; i < n - 2; i++) {
      // Load
      yi = y[i];
      yip2 = y[i+2];
      xi = x[i];
      xip1 = x[i+1];
      xip2 = x[i+2];
      zi = z[(i % 4) * 10]; // take outside of loop and calculate z[0], z[10], z[20], z[30]
      
      // Compute
      xi = xi * M_SQRT1_2 + yi * 0.2;
      xip1 = xip1 + zi * 0.3;
      xip2 = xip2 + sineValues[i % 6] * yip2; // take sin out and calculate the possible values 0, 60, 120 etc...
      
      // Store
      x[i]  = xi;
      x[i+1]  = xip1;
      x[i+2]  = xip2;
  }
}

void slow_performance7(float* x, float *y, float *z, int n) {
  float xi, xip1, xip2, xip3, xip4, xip5, xip6, xip7, xip8, xip9;
  float yi, yip1, yip2, yip3, yip4, yip5, yip6, yip7, yip8, yip9;
  float zi, zi1, zi2, zi3;
  float t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
  float t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36;
  float sin0, sin1, sin2, sin3, sin4, sin5;
  float c1 = 0.2;
  float c2 = 0.3;
  
  double mt_mpi_dt = 2 * M_PI / 3;
  float* sineValues = (float*) malloc(6 * sizeof(float));
  sineValues[0] = sin (0);
  sineValues[1] = sin (mt_mpi_dt);
  sineValues[2] = sin (mt_mpi_dt * 2);
  sineValues[3] = sin (mt_mpi_dt * 3);
  sineValues[4] = sin (mt_mpi_dt * 4);
  sineValues[5] = sin (mt_mpi_dt * 5);
  
  float m_sqrt1_2 = M_SQRT1_2;
  
  int limit = n - 7;
  int limit2 = n - 2;
  int i;
  for (i = 0; i < limit - 2; i+=8) {
      // Load 
      yi = y[i];
      yip1 = y[i+1];
      yip2 = y[i+2];
      yip3 = y[i+3];
      yip4 = y[i+4];
      yip5 = y[i+5];
      yip6 = y[i+6];
      yip7 = y[i+7];
      yip8 = y[i+8];
      yip9 = y[i+9];
      xi = x[i];
      xip1 = x[i+1];
      xip2 = x[i+2];
      xip3 = x[i+3];
      xip4 = x[i+4];
      xip5 = x[i+5];
      xip6 = x[i+6];
      xip7 = x[i+7];
      xip8 = x[i+8];
      xip9 = x[i+9];
      
      sin0 = sineValues[i%6];
      sin1 = sineValues[(i+1)%6];
      sin2 = sineValues[(i+2)%6];
      sin3 = sineValues[(i+3)%6];
      sin4 = sineValues[(i+4)%6];
      sin5 = sineValues[(i+5)%6];
      
      zi = z[(i%4) * 10] * c2;
      zi1 = z[((i+1)%4) * 10] * c2;
      zi2 = z[((i+2)%4) * 10] * c2;
      zi3 = z[((i+3)%4) * 10] * c2;
      
      // Compute
      t1 = yi * c1;
      t3 = sin0 * yip2;
      t4 = yip1 * c1;
      t6 = sin1 * yip3;
      t7 = yip2 * c1;
      t9 = sin2 * yip4;
      t10 = yip3 * c1;
      t12 = sin3 * yip5;
      t13 = yip4 * c1;
      t15 = sin4 * yip6;
      t16 = yip5 * c1;
      t18 = sin5 * yip7;
      t19 = yip6 * c1;
      t21 = sin0 * yip8;
      t22 = yip7 * c1;
      t24 = sin1 * yip9;
      
      xi = xi * m_sqrt1_2;
      xi = xi + t1;
      xip1 = xip1 + zi;
      xip2  = xip2 + t3;
      
      xip1 = xip1 * m_sqrt1_2;
      xip1 = xip1 + t4;
      xip2  = xip2 + zi1;
      xip3 = xip3 + t6;
      
      xip2 = xip2 * m_sqrt1_2;
      xip2 = xip2 + t7;
      xip3  = xip3 + zi2;
      xip4 = xip4 + t9;
      
      xip3 = xip3 * m_sqrt1_2;
      xip3  = xip3 + t10;
      xip4  = xip4 + zi3;
      xip5  = xip5 + t12;
      
      xip4 = xip4 * m_sqrt1_2;
      xip4 = xip4 + t13;
      xip5 = xip5 + zi;
      xip6 = xip6 + t15;
      
      xip5 = xip5 * m_sqrt1_2;
      xip5 = xip5 + t16;
      xip6 = xip6 + zi1;
      xip7 = xip7 + t18;
      
      xip6 = xip6 * m_sqrt1_2;
      xip6 = xip6 + t19;
      xip7 = xip7 + zi2;
      xip8 = xip8 + t21;
      
      xip7 = xip7 * m_sqrt1_2;
      xip7 = xip7 + t22;
      xip8 = xip8 + zi3;
      xip9 = xip9 + t24;
      
      // Store
      x[i] = xi;
      x[i+1] = xip1;
      x[i+2] = xip2;
      x[i+3] = xip3;
      x[i+4] = xip4;
      x[i+5] = xip5;
      x[i+6] = xip6;
      x[i+7] = xip7;
      x[i+8] = xip8;
      x[i+9] = xip9;
  }
  
  for (; i < n - 2; i++) {
      // Load
      yi = y[i];
      yip2 = y[i+2];
      xi = x[i];
      xip1 = x[i+1];
      xip2 = x[i+2];
      
      // Compute
      t1 = yi * c1;
      xi = xi * m_sqrt1_2;
      xi = xi + t1;
      zi = z[(i % 4) * 10] * c2;
      xip1 = xip1 + zi;
      t3 = sineValues[i % 6] * yip2;
      xip2 = xip2 + t3; // take sin out and calculate the possible values 0, 60, 120 etc...
      
      // Store
      x[i]  = xi;
      x[i+1]  = xip1;
      x[i+2]  = xip2;
  }
}

void slow_performance13(float* x, float *y, float *z, int n) {
  float xi, xip1, xip2, xip3, xip4, xip5, xip6, xip7, xip8, xip9, xip10, xip11, xip12, xip13;
  float yi, yip1, yip2, yip3, yip4, yip5, yip6, yip7, yip8, yip9, yip10, yip11, yip12, yip13;
  float zi, zi1, zi2, zi3;
  float t1;
  float sin1;
  float c1 = 0.2;
  float c2 = 0.3;
  zi = z[0] * c2;
  zi1 = z[10] * c2;
  zi2 = z[20] * c2;
  zi3 = z[30] * c2;
  
  sin1 = sin (M_PI / 3);
  float sineValues[3] = {0.0, sin1, -sin1};
  
  float m_sqrt1_2 = M_SQRT1_2;
  
  int i;
  for (i = 0; i < n - 13; i+=12) {
      // Load 
      yi = y[i];
      yip1 = y[i+1];
      yip2 = y[i+2];
      yip3 = y[i+3];
      yip4 = y[i+4];
      yip5 = y[i+5];
      yip6 = y[i+6];
      yip7 = y[i+7];
      yip8 = y[i+8];
      yip9 = y[i+9];
      yip10 = y[i+10];
      yip11 = y[i+11];
      yip12 = y[i+12];
      yip13 = y[i+13];
      xi = x[i];
      xip1 = x[i+1];
      xip2 = x[i+2];
      xip3 = x[i+3];
      xip4 = x[i+4];
      xip5 = x[i+5];
      xip6 = x[i+6];
      xip7 = x[i+7];
      xip8 = x[i+8];
      xip9 = x[i+9];
      xip10 = x[i+10];
      xip11 = x[i+11];
      xip12 = x[i+12];
      xip13 = x[i+13];
      
      // Compute
      xi = xi * m_sqrt1_2;
      t1 = yi * c1;
      x[i] = xi + t1;
      xip1 = xip1 + zi;
      
      xip1 = xip1 * m_sqrt1_2;
      t1 = yip1 * c1;
      x[i+1] = xip1 + t1;
      xip2  = xip2 + zi1;
      t1 = sin1 * yip3;
      xip3 = xip3 + t1;
      
      xip2 = xip2 * m_sqrt1_2;
      t1 = yip2 * c1;
      x[i+2] = xip2 + t1;
      xip3  = xip3 + zi2;
      t1 = sin1 * yip4;
      xip4 = xip4 - t1;
      
      xip3 = xip3 * m_sqrt1_2;
      t1 = yip3 * c1;
      x[i+3]  = xip3 + t1;
      xip4  = xip4 + zi3;
      
      xip4 = xip4 * m_sqrt1_2;
      t1 = yip4 * c1;
      x[i+4] = xip4 + t1;
      xip5 = xip5 + zi;
      t1 = sin1 * yip6;
      xip6 = xip6 + t1;
      
      xip5 = xip5 * m_sqrt1_2;
      t1 = yip5 * c1;
      x[i+5] = xip5 + t1;
      xip6 = xip6 + zi1;
      t1 = sin1 * yip7;
      xip7 = xip7 - t1;
      
      xip6 = xip6 * m_sqrt1_2;
      t1 = yip6 * c1;
      x[i+6] = xip6 + t1;
      xip7 = xip7 + zi2;
      
      xip7 = xip7 * m_sqrt1_2;
      t1 = yip7 * c1;
      x[i+7] = xip7 + t1;
      xip8 = xip8 + zi3;
      t1 = sin1 * yip9;
      xip9 = xip9 + t1;
      
      xip8 = xip8 * m_sqrt1_2;
      t1 = yip8 * c1;
      x[i+8] = xip8 + t1;
      xip9 = xip9 + zi;
      t1 = sin1 * yip10;
      xip10 = xip10 - t1;
      
      xip9 = xip9 * m_sqrt1_2;
      t1 = yip9 * c1;
      x[i+9] = xip9 + t1;
      xip10 = xip10 + zi1;
      
      xip10 = xip10 * m_sqrt1_2;
      t1 = yip10 * c1;
      x[i+10] = xip10 + t1;
      xip11 = xip11 + zi2;
      t1 = sin1 * yip12;
      xip12 = xip12 + t1;
      
      xip11 = xip11 * m_sqrt1_2;
      t1 = yip11 * c1;
      x[i+11] = xip11 + t1;
      x[i+12] = xip12 + zi3;
      t1 = sin1 * yip13;
      x[i+13] = xip13 - t1;
  }
  
  for (; i < n - 2; i++) {
      x[i]     = x[i] * m_sqrt1_2 + y[i] * c1;
      x[i+1]  += z[(i % 4) * 10] * c2;
      x[i+2]  += sineValues[i%3] * y[i+2];
  }
}

void maxperformance(float* x, float *y, float *z, int n) {
  float xip3, xip4, xip5, xip6, xip7, xip8, xip9, xip10;
  float yip3, yip4, yip6, yip7, yip9, yip10;
  float zi, zi1, zi2, zi3;
  float sin1;
  float c1 = 0.2;
  float c2 = 0.3;
  float m_sqrt1_2 = M_SQRT1_2;
  zi = z[0] * c2 * m_sqrt1_2;
  zi1 = z[10] * c2 * m_sqrt1_2;
  zi2 = z[20] * c2 * m_sqrt1_2;
  zi3 = z[30] * c2 * m_sqrt1_2;
  float zi3_2 = z[30] * c2;
  sin1 = sin (M_PI / 3);
  float sineValues[3] = {0.0, sin1, -sin1};
  
  float sin2 = sin1 * m_sqrt1_2 + c1;
  float sin3 = -sin1 * m_sqrt1_2 + c1;
  
  
  int i;
  for (i = 0; i < n - 13; i+=12) {
      // Load 
      yip3 = y[i+3];
      yip4 = y[i+4];
      yip6 = y[i+6];
      yip7 = y[i+7];
      yip9 = y[i+9];
      yip10 = y[i+10];
      xip3 = x[i+3];
      xip4 = x[i+4];
      xip5 = x[i+5];
      xip6 = x[i+6];
      xip7 = x[i+7];
      xip8 = x[i+8];
      xip9 = x[i+9];
      xip10 = x[i+10];
      
      // Compute
      x[i] = x[i] * m_sqrt1_2 + y[i] * c1;
      x[i+1] = x[i+1] * m_sqrt1_2 + y[i+1] * c1 + zi;
      x[i+2] = x[i+2] * m_sqrt1_2 + y[i+2] * c1 + zi1;
      x[i+3]  = xip3 * m_sqrt1_2 + yip3 * sin2 + zi2;
      x[i+4] = xip4 * m_sqrt1_2 + yip4 * sin3 + zi3;
      x[i+5] = xip5 * m_sqrt1_2 + y[i+5] * c1 + zi;
      x[i+6] = xip6 * m_sqrt1_2 + yip6 * sin2 + zi1;
      x[i+7] = xip7 * m_sqrt1_2 + yip7 * sin3 + zi2;
      x[i+8] = xip8 * m_sqrt1_2 + y[i+8] * c1 + zi3;
      x[i+9] = xip9 * m_sqrt1_2 + yip9 * sin2 + zi;
      x[i+10] = xip10 * m_sqrt1_2 + yip10 * sin3 + zi1;
      x[i+11] = x[i+11] * m_sqrt1_2 + y[i+11] * c1 + zi2;
      x[i+12] = x[i+12] + y[i+12] * sin1 + zi3_2;
      x[i+13] = x[i+13] - y[i+13] * sin1;
  }
  
  for (; i < n - 2; i++) {
      x[i]     = x[i] * m_sqrt1_2 + y[i] * c1;
      x[i+1]  += z[(i % 4) * 10] * c2;
      x[i+2]  += sineValues[i%3] * y[i+2];
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
  add_function(&slow_performance3, "slow_performance3",1);
  add_function(&slow_performance4, "slow_performance4",1);
  add_function(&slow_performance6, "slow_performance6",1);
  add_function(&slow_performance7, "slow_performance7",1);
  add_function(&slow_performance13, "slow_performance13",1);
  add_function(&maxperformance, "maxperformance",1);

}
