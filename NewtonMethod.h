#ifndef __NEWTONMETHOD__
#define __NEWTONMETHOD__

#include "polynomial.h"
#include <gmp.h>
#include <vector>
#include <iomanip>
class NewtonMethod
{
  private:
    mpf_class start_x;
    //for normal Newton
    polynomial u;
    polynomial derivate_u;

    // for Taylor series method
    mpf_class start_ux;
    mpf_class start_dux;
    vector<polynomial> coeff_ODE;
    vector<polynomial> d_coeff_ODE;
    vector<polynomial> d2_coeff_ODE;
  public:
    NewtonMethod() = default;
    //_initial_value = {u(start_x),u'(start_x)},coeff_ODE = {p,q,r}
    NewtonMethod(mpf_class _start_x, vector<mpf_class> _initial_value, vector<polynomial> _coeff_ODE);
    NewtonMethod(mpf_class _start_x, vector<polynomial> _initial_function, vector<polynomial> _coeff_ODE);
    //_initial_function = {u,u'}
    NewtonMethod(mpf_class _start_x, vector<polynomial> _initial_function);
    mpf_class compute_root_with_taylor(mpz_class taylor_order,mpz_class accuracy);
    mpf_class compute_root(mpz_class accuracy);
    void set_start_x(mpf_class _start_x);
    void set_start_value(vector<mpf_class> _initial_value);
    ~NewtonMethod() = default;
  private:
    mpf_class compute_eps(mpz_class accuracy);
    vector<mpf_class> compute_taylor_value(mpf_class current_x, mpf_class next_x, 
        vector<mpf_class> current_value,mpz_class tylor_order, unsigned int accuracy);
};

NewtonMethod::NewtonMethod(mpf_class _start_x, vector<mpf_class> _initial_value, vector<polynomial> _coeff_ODE)
  :start_x{_start_x},coeff_ODE{_coeff_ODE}
{
  start_ux = _initial_value.at(0);
  start_dux = _initial_value.at(1);
  int s = coeff_ODE.size();
  d_coeff_ODE.resize(s);
  d2_coeff_ODE.resize(s);
  for(int i = 0; i < s; ++i)
  {
    d_coeff_ODE.at(i) = derivate(coeff_ODE.at(i));
    d2_coeff_ODE.at(i) = derivate(d_coeff_ODE.at(i));
  }
};

NewtonMethod::NewtonMethod(mpf_class _start_x, vector<polynomial> _initial_function, vector<polynomial> _coeff_ODE)
  :start_x{_start_x},coeff_ODE{_coeff_ODE}
{
  u = _initial_function.at(0);
  derivate_u = _initial_function.at(1);
  start_ux = u(_start_x);
  start_dux = derivate_u(_start_x);
  int s = coeff_ODE.size();
  d_coeff_ODE.resize(s);
  d2_coeff_ODE.resize(s);
  for(int i = 0; i < s; ++i)
  {
    d_coeff_ODE.at(i) = derivate(coeff_ODE.at(i));
    d2_coeff_ODE.at(i) = derivate(d_coeff_ODE.at(i));
  }
};
NewtonMethod::NewtonMethod(mpf_class _start_x, vector<polynomial> _initial_function)
  :start_x{_start_x}
{
  u = _initial_function.at(0);
  derivate_u = _initial_function.at(1);
};

mpf_class NewtonMethod::compute_root_with_taylor(mpz_class taylor_order=30,mpz_class accuracy=16)
{
  if (coeff_ODE.size() == 0) {
    cout << "You don't have coeff_ODE vector, you should use compute_root function\n";
    return 0;
  }
  mpf_class eps = compute_eps(accuracy);
  cout << "eps is:\n";
  cout << fixed << setprecision(80) << eps  << endl;
  unsigned int acc = accuracy.get_ui();
  mpf_class current_x(start_x,2*acc);
  mpf_class current_ux(start_ux,2*acc);
  mpf_class current_dux(start_dux,2*acc);
  mpf_class next_x(0,2*acc);
  mpf_class h(0,2*acc);
  do {
    next_x = current_x - current_ux / current_dux;
    //compute u(x) and du(x) via taylor expansion.
    vector<mpf_class> current_value = {current_ux,current_dux};
    vector<mpf_class> next_value = compute_taylor_value(current_x, next_x, current_value, taylor_order, acc);
    current_ux = next_value.at(0);
    current_dux = next_value.at(1);
    current_x = next_x;
  } while(abs(u(current_x)) > eps); 
  //while(abs(current_ux) > eps);

  cout << "u(current_x) is:\n";
  cout << abs(u(current_x)) << endl;
  return current_x;
};

mpf_class NewtonMethod::compute_root(mpz_class accuracy=16)
{
  if (u.getDegree() == 0) {
    cout << "You don't have u function, you should use compute_root_with_taylor function\n";
    return 0;
  }
  mpf_class eps = compute_eps(accuracy); 
  cout << "eps is:\n";
  cout << fixed << setprecision(80) << eps << endl;
  unsigned int acc = accuracy.get_ui();
  mpf_class current_x(start_x,2*acc);
  mpf_class current_ux(u(current_x),2*acc);
  mpf_class current_dux(derivate_u(current_x),2*acc);
  mpf_class next_x(0,2*acc);

  do {
    next_x = current_x - current_ux / current_dux;
    current_x = next_x;
    current_ux = u(current_x);
    current_dux = derivate_u(current_x);
  } while( abs(current_ux) > eps);

  cout << "current_ux is:\n";
  cout << abs(current_ux) << endl;
  return current_x;
};
 
mpf_class NewtonMethod::compute_eps(mpz_class accuracy)
{
  unsigned int acc = accuracy.get_ui();
  mpf_t eps_t,base_t;
  mpf_init2(eps_t,2*acc);
  mpf_init2(base_t,2*acc);
  mpf_init_set_d(base_t,0.1);
  mpf_pow_ui(eps_t,base_t,acc);
  mpf_class eps(eps_t,2*acc);
  mpf_clear(eps_t);
  mpf_clear(base_t);
  return eps;
};

vector<mpf_class> NewtonMethod::compute_taylor_value(mpf_class current_x, mpf_class next_x, vector<mpf_class> current_value,mpz_class taylor_order, unsigned int accuracy)
{

  vector<mpf_class> current_ukx(taylor_order.get_ui()+1,mpf_class(0,accuracy));
  copy(cbegin(current_value),cend(current_value),begin(current_ukx));


  mpf_class px = coeff_ODE.at(0)(current_x);
  mpf_class qx = coeff_ODE.at(1)(current_x);
  mpf_class rx = coeff_ODE.at(2)(current_x);
  mpf_class p1x = d_coeff_ODE.at(0)(current_x);
  mpf_class q1x = d_coeff_ODE.at(1)(current_x);
  mpf_class r1x = d_coeff_ODE.at(2)(current_x);
  mpf_class p2x = d2_coeff_ODE.at(0)(current_x);
  mpf_class q2x = d2_coeff_ODE.at(1)(current_x);
  mpf_class r2x = d2_coeff_ODE.at(2)(current_x);
  for (int k = 0; k < current_ukx.size()-2; ++k)
  {
    //compute u(k+2)(x)
    mpf_class u_k2_x(0,2*accuracy);
    mpf_class coe1 = k*p1x + qx;
    if (coe1 != 0) {
      u_k2_x = u_k2_x - coe1*current_ukx.at(k+1);
    }
    mpf_class coe2 = 0.5*k*(k-1)*p2x + k*q1x + rx;
    if (coe2 != 0) {
      u_k2_x = u_k2_x - coe2*current_ukx.at(k);
    }
    mpf_class coe3 = 0.5*k*(k-1)*q2x + k*r1x;
    if (coe3 != 0) {
      u_k2_x = u_k2_x - coe3*current_ukx.at(k-1);
    }
    mpf_class coe4 = 0.5*k*(k-1)*r2x;
    if (coe4 != 0) {
      u_k2_x = u_k2_x - coe4*current_ukx.at(k-2);
    } 
    
    //BUG:if px == 0, what???
    current_ukx.at(k+2) = u_k2_x/px;
  }

  mpf_class next_ux(0,2*accuracy);
  mpf_class next_dux(0,2*accuracy);
  mpz_class factorial = 1;
  mpf_class hk = 1;
  mpf_class h(next_x - current_x,2*accuracy);
  for(int i = 0; i < taylor_order.get_ui(); ++i)
  {
    next_ux = next_ux + (current_ukx.at(i)/factorial)*hk;
    next_dux = next_dux + (current_ukx.at(i+1)/factorial)*hk;
    factorial = factorial*(i+1);
    hk = hk*h;
  }
  next_ux = next_ux + (current_ukx.back()/factorial)*hk;
  vector<mpf_class> next_value = {next_ux,next_dux};
  return next_value;
};
#else
//do nothing；
#endif
