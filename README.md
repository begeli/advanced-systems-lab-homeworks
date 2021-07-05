# Advanced Systems Lab Homeworks
My solutions for ETH Zurich Advanced Systems Lab homeworks. Homeworks were about analyzing performance of given code and applying optimization techniques learnt in class.

## The homework topics are as follows:
**Homework 1:** Cost, Performance and Operational Intensity Bound Anaylsis, ILP Analysis 

**Homework 2:** Applying standard C optimizations, Implementing Microbenchmarks to measure gap and latency of operations

* The optimized code managed to speedup the _unoptimized_ implementation by **20.21** times

**Homework 3:** Vectorization with Intel AVX intrinsics (Optimizing FIR filters, vector addition and complex conversion)

* The vectorized implementation of FIR filters managed to speedup the _optimized_ implementation by **2.666** times
* The vectorized implementation of vector addition managed to speedup the _optimized_ implementation by **3.568** times

**Homework 4:** Cache Associativity, Cachme Miss Analysis, Roofline Models

## Homework platform

All the code was optimized for the `Intel Xeon Silver 4210 Processor`. The platform I tested the code on compiled the code using `GCC 8.3.1` and flags 
`-O3 -fno-tree-vectorize -march=skylake`.
