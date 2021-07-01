#ifndef MICROBENCHMARK_H
#define MICROBENCHMARK_H

typedef enum {
    ADD_LAT,
    ADD_GAP,
    DIV_LAT,
    DIV_GAP,
    DIV_LAT_MIN,
    DIV_GAP_MIN,
    FOO_LAT,
    FOO_GAP,
    FOO_LAT_MIN,
    FOO_GAP_MIN,
    
    /* Limits */
    START_TEST = ADD_LAT,
    END_TEST = FOO_GAP_MIN,
} microbenchmark_mode_t;


void    initialize_microbenchmark_data(microbenchmark_mode_t mode);
double  microbenchmark_get_add_latency ();
double  microbenchmark_get_add_gap     ();
double  microbenchmark_get_div_latency ();
double  microbenchmark_get_div_gap     ();
double  microbenchmark_get_sqrt_latency();
double  microbenchmark_get_sqrt_gap    ();
double  microbenchmark_get_foo_latency();
double  microbenchmark_get_foo_gap    ();


#endif /* MICROBENCHMARK_H */
