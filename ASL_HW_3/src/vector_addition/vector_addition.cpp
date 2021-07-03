#include "common.h"
#include <immintrin.h>
#include <math.h>

void slow_performance1(vector_t x, vector_t y, vector_t z, int n) {
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

void slow_performance3(vector_t x, vector_t y, vector_t z, int n) {
  int i;
  __m256d xValv, yValv;
  __m256d mxValv, myValv;
  __m256d xAnd, yAndNot, xFabAdd, yFabAdd, newAddVals;
  __m256d idNotEqRes, idEqRes;
  __m256i xIdv, yIdv, zIdv;
  
  __m256d idCmp, fabCmp;
  
  __m256d totSum = _mm256_set1_pd(0.0);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  for (i = 0; i < n - 3; i+=4) {
    // compare for their equality, store in new vector
    xIdv = _mm256_load_si256((__m256i *) (x.id + i));
    yIdv = _mm256_load_si256((__m256i *) (y.id + i));
    xValv = _mm256_load_pd(x.val + i);
    yValv = _mm256_load_pd(y.val + i); 
    idEqRes = _mm256_add_pd(xValv, yValv); 
    
    idCmp = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv, yIdv));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv = _mm256_and_pd(xValv, maskAbs); 
    myValv = _mm256_and_pd(yValv, maskAbs); 
    
    fabCmp = _mm256_cmp_pd(mxValv, myValv, _CMP_GT_OQ); 
    
    // If fabs(x.val[i]) > fabs(y.val[i])
    xAnd = _mm256_and_pd(xValv, fabCmp); 
    yFabAdd = _mm256_and_pd(myValv, fabCmp); 
    xIdv = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv), fabCmp)); 
    
    // If fabs(x.val[i]) <= fabs(y.val[i])
    yAndNot = _mm256_andnot_pd(fabCmp, yValv); 
    xFabAdd = _mm256_andnot_pd(fabCmp, mxValv); 
    yIdv = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp, _mm256_castsi256_pd(yIdv))); 
    
    // Combine Everything
    idNotEqRes = _mm256_or_pd(xAnd, yAndNot);
    newAddVals = _mm256_or_pd(xFabAdd, yFabAdd);
    zIdv = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv), _mm256_castsi256_pd(yIdv)));
    
    idNotEqRes = _mm256_andnot_pd(idCmp, idNotEqRes);
    idEqRes = _mm256_and_pd(idCmp, idEqRes);
    newAddVals = _mm256_andnot_pd(idCmp, newAddVals);
    totSum = _mm256_add_pd(totSum, newAddVals);
    
    // Store results
    _mm256_store_pd(z.val + i, _mm256_or_pd(idEqRes, idNotEqRes));
    _mm256_store_si256((__m256i *)(z.id + i), zIdv);
  }
  
  double sum[4];
  _mm256_storeu_pd(sum, totSum);
  z.val[n] += sum[0] + sum[1] + sum[2] + sum[3];
  
  for (; i < n; i++) {
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

void slow_performance4(vector_t x, vector_t y, vector_t z, int n) {
  int i;
  __m256d xValv, yValv;
  __m256d mxValv, myValv;
  __m256d xAnd, yAndNot, xFabAdd, yFabAdd, newAddVals;
  __m256d idNotEqRes, idEqRes;
  __m256i xIdv, yIdv, zIdv;
  
  __m256d idCmp, fabCmp;
  
  __m256d totSum = _mm256_set1_pd(0.0);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  __m256d xValv2, yValv2;
  __m256d mxValv2, myValv2;
  __m256d xAnd2, yAndNot2, xFabAdd2, yFabAdd2, newAddVals2;
  __m256d idNotEqRes2, idEqRes2;
  __m256i xIdv2, yIdv2, zIdv2;
  
  __m256d idCmp2, fabCmp2;
  
  __m256d totSum2 = _mm256_set1_pd(0.0);
  __m256d maskAbs2 = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  __m256d xValv3, yValv3;
  __m256d mxValv3, myValv3;
  __m256d xAnd3, yAndNot3, xFabAdd3, yFabAdd3, newAddVals3;
  __m256d idNotEqRes3, idEqRes3;
  __m256i xIdv3, yIdv3, zIdv3;
  
  __m256d idCmp3, fabCmp3;
  
  __m256d totSum3 = _mm256_set1_pd(0.0);
  __m256d maskAbs3 = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  __m256d xValv4, yValv4;
  __m256d mxValv4, myValv4;
  __m256d xAnd4, yAndNot4, xFabAdd4, yFabAdd4, newAddVals4;
  __m256d idNotEqRes4, idEqRes4;
  __m256i xIdv4, yIdv4, zIdv4;
  
  __m256d idCmp4, fabCmp4;
  
  __m256d totSum4 = _mm256_set1_pd(0.0);
  __m256d maskAbs4 = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  for (i = 0; i < n - 15; i+=16) {
    // compare for their equality, store in new vector
    xIdv = _mm256_load_si256((__m256i *) (x.id + i));
    yIdv = _mm256_load_si256((__m256i *) (y.id + i));
    xValv = _mm256_load_pd(x.val + i);
    yValv = _mm256_load_pd(y.val + i); 
    idEqRes = _mm256_add_pd(xValv, yValv); 
    
    idCmp = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv, yIdv));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv = _mm256_and_pd(xValv, maskAbs); 
    myValv = _mm256_and_pd(yValv, maskAbs); 
    
    fabCmp = _mm256_cmp_pd(mxValv, myValv, _CMP_GT_OQ); 
    
    // If fabs(x.val[i]) > fabs(y.val[i])
    xAnd = _mm256_and_pd(xValv, fabCmp); 
    yFabAdd = _mm256_and_pd(myValv, fabCmp); 
    xIdv = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv), fabCmp)); 
    
    // If fabs(x.val[i]) <= fabs(y.val[i])
    yAndNot = _mm256_andnot_pd(fabCmp, yValv); 
    xFabAdd = _mm256_andnot_pd(fabCmp, mxValv); 
    yIdv = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp, _mm256_castsi256_pd(yIdv))); 
    
    // Combine Everything
    idNotEqRes = _mm256_or_pd(xAnd, yAndNot);
    newAddVals = _mm256_or_pd(xFabAdd, yFabAdd);
    zIdv = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv), _mm256_castsi256_pd(yIdv)));
    
    idNotEqRes = _mm256_andnot_pd(idCmp, idNotEqRes);
    idEqRes = _mm256_and_pd(idCmp, idEqRes);
    newAddVals = _mm256_andnot_pd(idCmp, newAddVals);
    totSum = _mm256_add_pd(totSum, newAddVals);
    
    // Store results
    _mm256_store_pd(z.val + i, _mm256_or_pd(idEqRes, idNotEqRes));
    _mm256_store_si256((__m256i *)(z.id + i), zIdv);
    
    // Unroll 1
    xIdv2 = _mm256_load_si256((__m256i *) (x.id + i + 4));
    yIdv2 = _mm256_load_si256((__m256i *) (y.id + i + 4));
    xValv2 = _mm256_load_pd(x.val + i + 4);
    yValv2 = _mm256_load_pd(y.val + i + 4); 
    idEqRes2 = _mm256_add_pd(xValv2, yValv2); 
    
    idCmp2 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv2, yIdv2));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv2 = _mm256_and_pd(xValv2, maskAbs2); 
    myValv2 = _mm256_and_pd(yValv2, maskAbs2); 
    
    fabCmp2 = _mm256_cmp_pd(mxValv2, myValv2, _CMP_GT_OQ); 
    
    // If fabs(x.val[i]) > fabs(y.val[i])
    xAnd2 = _mm256_and_pd(xValv2, fabCmp2); 
    yFabAdd2 = _mm256_and_pd(myValv2, fabCmp2); 
    xIdv2 = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv2), fabCmp2)); 
    
    // If fabs(x.val[i]) <= fabs(y.val[i])
    yAndNot2 = _mm256_andnot_pd(fabCmp2, yValv2); 
    xFabAdd2 = _mm256_andnot_pd(fabCmp2, mxValv2); 
    yIdv2 = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp2, _mm256_castsi256_pd(yIdv2))); 
    
    // Combine Everything
    idNotEqRes2 = _mm256_or_pd(xAnd2, yAndNot2);
    newAddVals2 = _mm256_or_pd(xFabAdd2, yFabAdd2);
    zIdv2 = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv2), _mm256_castsi256_pd(yIdv2)));
    
    idNotEqRes2 = _mm256_andnot_pd(idCmp2, idNotEqRes2);
    idEqRes2 = _mm256_and_pd(idCmp2, idEqRes2);
    newAddVals2 = _mm256_andnot_pd(idCmp2, newAddVals2);
    totSum2 = _mm256_add_pd(totSum2, newAddVals2);
    
    // Store results
    _mm256_store_pd(z.val + i + 4, _mm256_or_pd(idEqRes2, idNotEqRes2));
    _mm256_store_si256((__m256i *)(z.id + i + 4), zIdv2);
    
    // Unroll 2
    xIdv3 = _mm256_load_si256((__m256i *) (x.id + i + 8));
    yIdv3 = _mm256_load_si256((__m256i *) (y.id + i + 8));
    xValv3 = _mm256_load_pd(x.val + i + 8);
    yValv3 = _mm256_load_pd(y.val + i + 8); 
    idEqRes3 = _mm256_add_pd(xValv3, yValv3); 
    
    idCmp3 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv3, yIdv3));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv3 = _mm256_and_pd(xValv3, maskAbs3); 
    myValv3 = _mm256_and_pd(yValv3, maskAbs3); 
    
    fabCmp3 = _mm256_cmp_pd(mxValv3, myValv3, _CMP_GT_OQ); 
    
    // If fabs(x.val[i]) > fabs(y.val[i])
    xAnd3 = _mm256_and_pd(xValv3, fabCmp3); 
    yFabAdd3 = _mm256_and_pd(myValv3, fabCmp3); 
    xIdv3 = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv3), fabCmp3)); 
    
    // If fabs(x.val[i]) <= fabs(y.val[i])
    yAndNot3 = _mm256_andnot_pd(fabCmp3, yValv3); 
    xFabAdd3 = _mm256_andnot_pd(fabCmp3, mxValv3); 
    yIdv3 = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp3, _mm256_castsi256_pd(yIdv3))); 
    
    // Combine Everything
    idNotEqRes3 = _mm256_or_pd(xAnd3, yAndNot3);
    newAddVals3 = _mm256_or_pd(xFabAdd3, yFabAdd3);
    zIdv3 = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv3), _mm256_castsi256_pd(yIdv3)));
    
    idNotEqRes3 = _mm256_andnot_pd(idCmp3, idNotEqRes3);
    idEqRes3 = _mm256_and_pd(idCmp3, idEqRes3);
    newAddVals3 = _mm256_andnot_pd(idCmp3, newAddVals3);
    totSum3 = _mm256_add_pd(totSum3, newAddVals3);
    
    // Store results
    _mm256_store_pd(z.val + i + 8, _mm256_or_pd(idEqRes3, idNotEqRes3));
    _mm256_store_si256((__m256i *)(z.id + i + 8), zIdv3);
    
    // Unroll 3
    xIdv4 = _mm256_load_si256((__m256i *) (x.id + i + 12));
    yIdv4 = _mm256_load_si256((__m256i *) (y.id + i + 12));
    xValv4 = _mm256_load_pd(x.val + i + 12);
    yValv4 = _mm256_load_pd(y.val + i + 12); 
    idEqRes4 = _mm256_add_pd(xValv4, yValv4); 
    
    idCmp4 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv4, yIdv4));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv4 = _mm256_and_pd(xValv4, maskAbs4); 
    myValv4 = _mm256_and_pd(yValv4, maskAbs4); 
    
    fabCmp4 = _mm256_cmp_pd(mxValv4, myValv4, _CMP_GT_OQ); 
    
    // If fabs(x.val[i]) > fabs(y.val[i])
    xAnd4 = _mm256_and_pd(xValv4, fabCmp4); 
    yFabAdd4 = _mm256_and_pd(myValv4, fabCmp4); 
    xIdv4 = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv4), fabCmp4)); 
    
    // If fabs(x.val[i]) <= fabs(y.val[i])
    yAndNot4 = _mm256_andnot_pd(fabCmp4, yValv4); 
    xFabAdd4 = _mm256_andnot_pd(fabCmp4, mxValv4); 
    yIdv4 = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp4, _mm256_castsi256_pd(yIdv4))); 
    
    // Combine Everything
    idNotEqRes4 = _mm256_or_pd(xAnd4, yAndNot4);
    newAddVals4 = _mm256_or_pd(xFabAdd4, yFabAdd4);
    zIdv4 = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv4), _mm256_castsi256_pd(yIdv4)));
    
    idNotEqRes4 = _mm256_andnot_pd(idCmp4, idNotEqRes4);
    idEqRes4 = _mm256_and_pd(idCmp4, idEqRes4);
    newAddVals4 = _mm256_andnot_pd(idCmp4, newAddVals4);
    totSum4 = _mm256_add_pd(totSum4, newAddVals4);
    
    // Store results
    _mm256_store_pd(z.val + i + 12, _mm256_or_pd(idEqRes4, idNotEqRes4));
    _mm256_store_si256((__m256i *)(z.id + i + 12), zIdv4);
  }
  
  double sum[4];
  _mm256_storeu_pd(sum, totSum);
  z.val[n] += sum[0] + sum[1] + sum[2] + sum[3];
  _mm256_storeu_pd(sum, totSum2);
  z.val[n] += sum[0] + sum[1] + sum[2] + sum[3];
  _mm256_storeu_pd(sum, totSum3);
  z.val[n] += sum[0] + sum[1] + sum[2] + sum[3];
  _mm256_storeu_pd(sum, totSum4);
  z.val[n] += sum[0] + sum[1] + sum[2] + sum[3];

  for (; i < n; i++) {
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

void maxperformance(vector_t x, vector_t y, vector_t z, int n) {
  int i;
  __m256d xValv, yValv;
  __m256d mxValv, myValv;
  __m256d xAnd, yAndNot, xFabAdd, yFabAdd, newAddVals;
  __m256d idNotEqRes, idEqRes;
  __m256i xIdv, yIdv, zIdv;
  
  __m256d idCmp, fabCmp;
  
  __m256d totSum = _mm256_set1_pd(0.0);
  __m256d maskAbs = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  __m256d xValv2, yValv2;
  __m256d mxValv2, myValv2;
  __m256d xAnd2, yAndNot2, xFabAdd2, yFabAdd2, newAddVals2;
  __m256d idNotEqRes2, idEqRes2;
  __m256i xIdv2, yIdv2, zIdv2;
  
  __m256d idCmp2, fabCmp2;
  
  __m256d totSum2 = _mm256_set1_pd(0.0);
  __m256d maskAbs2 = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  __m256d xValv3, yValv3;
  __m256d mxValv3, myValv3;
  __m256d xAnd3, yAndNot3, xFabAdd3, yFabAdd3, newAddVals3;
  __m256d idNotEqRes3, idEqRes3;
  __m256i xIdv3, yIdv3, zIdv3;
  
  __m256d idCmp3, fabCmp3;
  
  __m256d totSum3 = _mm256_set1_pd(0.0);
  __m256d maskAbs3 = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  __m256d xValv4, yValv4;
  __m256d mxValv4, myValv4;
  __m256d xAnd4, yAndNot4, xFabAdd4, yFabAdd4, newAddVals4;
  __m256d idNotEqRes4, idEqRes4;
  __m256i xIdv4, yIdv4, zIdv4;
  
  __m256d idCmp4, fabCmp4;
  
  __m256d totSum4 = _mm256_set1_pd(0.0);
  __m256d maskAbs4 = _mm256_castsi256_pd(_mm256_set1_epi64x(0x7FFFFFFFFFFFFFFF));
  
  for (i = 0; i < n - 15; i+=16) {
    // compare for their equality, store in new vector
    xIdv = _mm256_load_si256((__m256i *) (x.id + i));
    yIdv = _mm256_load_si256((__m256i *) (y.id + i));
    xValv = _mm256_load_pd(x.val + i);
    yValv = _mm256_load_pd(y.val + i); 
    idEqRes = _mm256_add_pd(xValv, yValv); 
    
    idCmp = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv, yIdv));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv = _mm256_and_pd(xValv, maskAbs); 
    myValv = _mm256_and_pd(yValv, maskAbs); 
    
    fabCmp = _mm256_cmp_pd(mxValv, myValv, _CMP_GT_OQ); 
    
    xIdv = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv), fabCmp)); 
    yIdv = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp, _mm256_castsi256_pd(yIdv))); 
    
    // Combine Everything
    idNotEqRes = _mm256_blendv_pd(yValv, xValv, fabCmp);
    newAddVals = _mm256_blendv_pd(mxValv, myValv, fabCmp);
    zIdv = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv), _mm256_castsi256_pd(yIdv)));
    
    newAddVals = _mm256_andnot_pd(idCmp, newAddVals);
    totSum = _mm256_add_pd(totSum, newAddVals);
    
    _mm256_store_pd(z.val + i, _mm256_blendv_pd(idNotEqRes, idEqRes, idCmp));
    _mm256_store_si256((__m256i *)(z.id + i), zIdv);
    
    // Unroll 1
    xIdv2 = _mm256_load_si256((__m256i *) (x.id + i + 4));
    yIdv2 = _mm256_load_si256((__m256i *) (y.id + i + 4));
    xValv2 = _mm256_load_pd(x.val + i + 4);
    yValv2 = _mm256_load_pd(y.val + i + 4); 
    idEqRes2 = _mm256_add_pd(xValv2, yValv2); 
    
    idCmp2 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv2, yIdv2));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv2 = _mm256_and_pd(xValv2, maskAbs2); 
    myValv2 = _mm256_and_pd(yValv2, maskAbs2); 
    
    fabCmp2 = _mm256_cmp_pd(mxValv2, myValv2, _CMP_GT_OQ); 
    
    xIdv2 = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv2), fabCmp2)); 
    yIdv2 = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp2, _mm256_castsi256_pd(yIdv2))); 
    
    // Combine Everything
    idNotEqRes2 = _mm256_blendv_pd(yValv2, xValv2, fabCmp2);
    newAddVals2 = _mm256_blendv_pd(mxValv2, myValv2, fabCmp2);
    zIdv2 = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv2), _mm256_castsi256_pd(yIdv2)));
    
    newAddVals2 = _mm256_andnot_pd(idCmp2, newAddVals2);
    totSum2 = _mm256_add_pd(totSum2, newAddVals2);
    
    // Store results
    _mm256_store_pd(z.val + i + 4, _mm256_blendv_pd(idNotEqRes2, idEqRes2, idCmp2));
    _mm256_store_si256((__m256i *)(z.id + i + 4), zIdv2);
    
    // Unroll 2
    xIdv3 = _mm256_load_si256((__m256i *) (x.id + i + 8));
    yIdv3 = _mm256_load_si256((__m256i *) (y.id + i + 8));
    xValv3 = _mm256_load_pd(x.val + i + 8);
    yValv3 = _mm256_load_pd(y.val + i + 8); 
    idEqRes3 = _mm256_add_pd(xValv3, yValv3); 
    
    idCmp3 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv3, yIdv3));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv3 = _mm256_and_pd(xValv3, maskAbs3); 
    myValv3 = _mm256_and_pd(yValv3, maskAbs3); 
    
    fabCmp3 = _mm256_cmp_pd(mxValv3, myValv3, _CMP_GT_OQ); 
    
    xIdv3 = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv3), fabCmp3)); 
    yIdv3 = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp3, _mm256_castsi256_pd(yIdv3))); 
    
    // Combine Everything
    idNotEqRes3 = _mm256_blendv_pd(yValv3, xValv3, fabCmp3);
    newAddVals3 = _mm256_blendv_pd(mxValv3, myValv3, fabCmp3);
    zIdv3 = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv3), _mm256_castsi256_pd(yIdv3)));
  
    newAddVals3 = _mm256_andnot_pd(idCmp3, newAddVals3);
    totSum3 = _mm256_add_pd(totSum3, newAddVals3);
    
    // Store results
    _mm256_store_pd(z.val + i + 8, _mm256_blendv_pd(idNotEqRes3, idEqRes3, idCmp3));
    _mm256_store_si256((__m256i *)(z.id + i + 8), zIdv3);
    
    // Unroll 3
    xIdv4 = _mm256_load_si256((__m256i *) (x.id + i + 12));
    yIdv4 = _mm256_load_si256((__m256i *) (y.id + i + 12));
    xValv4 = _mm256_load_pd(x.val + i + 12);
    yValv4 = _mm256_load_pd(y.val + i + 12); 
    idEqRes4 = _mm256_add_pd(xValv4, yValv4); 
    
    idCmp4 = _mm256_castsi256_pd(_mm256_cmpeq_epi64(xIdv4, yIdv4));
    
    // If Ids aren't equal
    /*---------------------------------------------------*/
    mxValv4 = _mm256_and_pd(xValv4, maskAbs4); 
    myValv4 = _mm256_and_pd(yValv4, maskAbs4); 
    
    fabCmp4 = _mm256_cmp_pd(mxValv4, myValv4, _CMP_GT_OQ); 
    
    // If fabs(x.val[i]) > fabs(y.val[i])
    xAnd4 = _mm256_and_pd(xValv4, fabCmp4); 
    yFabAdd4 = _mm256_and_pd(myValv4, fabCmp4); 
    xIdv4 = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(xIdv4), fabCmp4)); 
    
    // If fabs(x.val[i]) <= fabs(y.val[i])
    yAndNot4 = _mm256_andnot_pd(fabCmp4, yValv4); 
    xFabAdd4 = _mm256_andnot_pd(fabCmp4, mxValv4); 
    yIdv4 = _mm256_castpd_si256(_mm256_andnot_pd(fabCmp4, _mm256_castsi256_pd(yIdv4)));
    
    // Combine Everything
    idNotEqRes4 = _mm256_or_pd(xAnd4, yAndNot4);
    newAddVals4 = _mm256_or_pd(xFabAdd4, yFabAdd4);
    zIdv4 = _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(xIdv4), _mm256_castsi256_pd(yIdv4)));
    
    idNotEqRes4 = _mm256_andnot_pd(idCmp4, idNotEqRes4);
    idEqRes4 = _mm256_and_pd(idCmp4, idEqRes4);
    newAddVals4 = _mm256_andnot_pd(idCmp4, newAddVals4);
    totSum4 = _mm256_add_pd(totSum4, newAddVals4);
    
    // Store results
    _mm256_store_pd(z.val + i + 12, _mm256_or_pd(idEqRes4, idNotEqRes4));
    _mm256_store_si256((__m256i *)(z.id + i + 12), zIdv4);
  }
  
  double sum[4];
  totSum = _mm256_add_pd(totSum, totSum2);
  totSum = _mm256_add_pd(totSum, totSum3);
  totSum = _mm256_add_pd(totSum, totSum4);
  _mm256_storeu_pd(sum, totSum);
  z.val[n] += sum[0] + sum[1] + sum[2] + sum[3];
  
  for (; i < n; i++) {
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

/*
* Called by the driver to register your functions
* Use add_function(func, description) to add your own functions
*/
void register_functions() 
{
  add_function(&slow_performance1, "slow_performance1",1);
  add_function(&slow_performance3, "slow_performance3",1);
  add_function(&slow_performance4, "slow_performance4",1);
  add_function(&maxperformance, "maxperformance",1);
}
