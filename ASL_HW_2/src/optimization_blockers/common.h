#include <string>
typedef void(*comp_func)(float *, float *, float*, int);
void add_function(comp_func f, std::string name, int flop);
void kernel_base(float* A, float *x, float *y, int);