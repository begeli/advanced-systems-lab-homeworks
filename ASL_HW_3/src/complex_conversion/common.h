#include <string>
typedef void(*comp_func)(double *, double *, double*, int);
void add_function(comp_func f, std::string name, int flop);
void kernel_base(double*, double *, double *, int);