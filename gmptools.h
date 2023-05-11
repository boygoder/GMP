#ifndef __GMPTOOLS__
#define __GMPTOOLS__

#include <gmpxx.h>

mpf_class compute_pi(mpz_class order);

mpf_class mpf_class_pow_ui(mpf_class base,unsigned int exception);

const unsigned int precision = 512;
//低精度，大概十几位。
//const mpf_class pi(4.0*atan(1.0),precision);
const mpf_class pi  = compute_pi(256);
// if num<1e-160, we think num is zero.
const mpf_class zero_eps = mpf_class_pow_ui(mpf_class(0.1,precision), 160);
// if num > 1e160, we think num is infinity.
const mpf_class infinity = mpf_class_pow_ui(mpf_class(10,precision), 160);

#else
//do nothing
#endif
