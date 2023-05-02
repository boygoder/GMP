#ifndef __GMPTOOLS__
#define __GMPTOOLS__

#include <gmpxx.h>
#include <cmath>
const unsigned int precision = 256;
const mpf_class pi(4.0*atan(1.0),precision);




mpf_class mpf_class_pow_ui(mpf_class base,unsigned int exception)
{
  mpf_t base_t,result_t;
  mpf_init2(base_t,precision);
  mpf_set(base_t,base.get_mpf_t());
  mpf_init2(result_t,precision);
  mpf_pow_ui(result_t,base_t,exception);
  mpf_class result(result_t,precision);
  mpf_clear(base_t);
  mpf_clear(result_t);
  return result;
};

// if num<1e-80, we think num is zero.
const mpf_class zero_eps = mpf_class_pow_ui(mpf_class(0.1,precision), 80);
const mpf_class infinity = mpf_class_pow_ui(mpf_class(10,256), 80);
#else
//do nothing
#endif
