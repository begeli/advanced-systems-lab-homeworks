double scalar_product(double* x, double *y, int n) {
  /* Insert here your optimized implementation. */
    /*
      Latency od addition and multiplication operations in my computer is 4 cycles.
      Throughput of the addtion and multiplication operations 
      in my computer is 2 ops/cycle.
    */
    const int OPTIMAL_UNROLL_AMOUNT = 8; 
    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    double s3 = 0.0;
    double s4 = 0.0;
    double s5 = 0.0;
    double s6 = 0.0;
    double s7 = 0.0;

    int i;
    int limit = n - OPTIMAL_UNROLL_AMOUNT + 1;

    for (i = 0; i < limit; i += OPTIMAL_UNROLL_AMOUNT) {
        s0 += x[i] * y[i];
        s1 += x[i + 1] * y[i + 1];
        s2 += x[i + 2] * y[i + 2];
        s3 += x[i + 3] * y[i + 3];
        s4 += x[i + 4] * y[i + 4];
        s5 += x[i + 5] * y[i + 5];
        s6 += x[i + 6] * y[i + 6];
        s7 += x[i + 7] * y[i + 7];
    }

    // shouldn't happen since size is a multiple of 2^3
    for (; i < n; i++) {
        s0 += x[i] * y[i];
    }

  return s0 + s1 + s2 + s3 + s4 + s5 + s6 + s7;
}