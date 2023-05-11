#include "gmptools.h"
#include <cmath>
#include <vector>
#include <numeric>
using namespace std;
mpf_class compute_pi(mpz_class order)
{
  vector<mpf_class> Rn(2*order.get_ui(),mpf_class(0,precision));
  Rn.at(1) = 1;
  for(int i = 2; i < Rn.size(); ++i) {
    mpf_class tmp = mpf_class(i-1,precision)/mpf_class(2*i-1,precision); 
    Rn.at(i) = Rn.at(i-1) * tmp;
  }
  mpf_class half_pi = accumulate(Rn.begin(),Rn.end(),mpf_class(0,precision));
  mpf_class pi = 2.0 * half_pi;
  return pi;
};

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
