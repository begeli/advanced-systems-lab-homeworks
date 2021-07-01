#include <string>

struct vector_t {
  double*   val;
  uint64_t*  id;
};

typedef void(*comp_func)(vector_t, vector_t, vector_t, int);
void add_function(comp_func f, std::string name, int flop);
void kernel_base(vector_t, vector_t, vector_t, int);
