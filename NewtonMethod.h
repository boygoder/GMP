#ifndef __NEWTONMETHOD__
#define __NEWTONMETHOD__

#include "polynomial.h"
#include <vector>

class NewtonMethod
{
  private:
    mpf_class start_x{0,precision};
    //for normal Newton
    polynomial u;
    polynomial derivate_u;

    // for Taylor series method
    vector<polynomial> coeff_ODE;
    vector<polynomial> d_coeff_ODE;
    vector<polynomial> d2_coeff_ODE;
  public:
    NewtonMethod() = default;
    NewtonMethod(mpf_class _start_x, vector<polynomial> _initial_function, vector<polynomial> _coeff_ODE);
    //_initial_function = {u,u'}
    NewtonMethod(mpf_class _start_x, vector<polynomial> _initial_function);
    mpf_class compute_root_with_taylor(mpz_class taylor_order,mpz_class accuracy);
    mpf_class compute_root(mpz_class accuracy);
    void set_start_x(mpf_class _start_x);
    ~NewtonMethod() = default;
  private:
    mpf_class compute_eps(mpz_class accuracy);
    vector<mpf_class> compute_derivates(mpf_class current_x,mpz_class order);
    vector<mpf_class> compute_taylor_value(mpf_class current_x, mpf_class next_x,vector<mpf_class>& current_ukx,mpz_class tylor_order);
};




#else
//do nothing；
#endif
