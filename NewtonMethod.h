#ifndef __NEWTONMETHOD__
#define __NEWTONMETHOD__

#include "polynomial.h"
#include <vector>
#include <iomanip>
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



void NewtonMethod::set_start_x(mpf_class _start_x)
{
  start_x = _start_x;
};

NewtonMethod::NewtonMethod(mpf_class _start_x, vector<polynomial> _initial_function, vector<polynomial> _coeff_ODE)
  :start_x{_start_x},coeff_ODE{_coeff_ODE}
{
  u = _initial_function.at(0);
  derivate_u = _initial_function.at(1);
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
  mpf_class current_x(start_x,precision);
  mpf_class current_ux = u(current_x);
  mpf_class current_dux = derivate_u(current_x);
  vector<mpf_class> start_ukx = compute_derivates(current_x, taylor_order);
  mpf_class next_x(0,precision);
  mpf_class h(0,precision);
  do {
    next_x = current_x - current_ux / current_dux;
    //compute u(x) and du(x) via taylor expansion.
    //change current_x to start_x
    vector<mpf_class> next_value = compute_taylor_value(start_x, next_x,start_ukx,taylor_order);
    current_ux = next_value.at(0);
    current_dux = next_value.at(1);
    current_x = next_x;
    //cout << "u(current_x) is:\n";
    //cout << abs(u(current_x)) << endl;
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
  mpf_class current_x(start_x,precision);
  mpf_class current_ux(u(current_x),precision);
  mpf_class current_dux(derivate_u(current_x),precision);
  mpf_class next_x(0,precision);

  do {
    next_x = current_x - current_ux / current_dux;
    current_x = next_x;
    current_ux = u(current_x);
    current_dux = derivate_u(current_x);

    //cout << "current_ux is:\n";
    //cout << abs(current_ux) << endl;
  } while( abs(current_ux) > eps);

  cout << "current_ux is:\n";
  cout << abs(current_ux) << endl;
  return current_x;
};

mpf_class NewtonMethod::compute_eps(mpz_class accuracy)
{
  unsigned int acc = accuracy.get_ui();
  mpf_t eps_t,base_t;
  mpf_init2(eps_t,precision);
  mpf_init2(base_t,precision);
  mpf_class base(0.1,precision);
  mpf_set(base_t,base.get_mpf_t());
  mpf_pow_ui(eps_t,base_t,acc);
  mpf_class eps(eps_t,precision);
  mpf_clear(eps_t);
  mpf_clear(base_t);
  return eps;
};
vector<mpf_class> NewtonMethod::compute_derivates(mpf_class current_x,mpz_class order)
{
  vector<mpf_class> current_ukx(order.get_ui()+1,mpf_class(0,precision));
  current_ukx.at(0) = u(current_x);
  current_ukx.at(1) = derivate_u(current_x);

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
    mpf_class u_k2_x(0,precision);
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
  return current_ukx;

};
vector<mpf_class> NewtonMethod::compute_taylor_value(mpf_class current_x, mpf_class next_x,vector<mpf_class>& start_ukx,mpz_class taylor_order)
{
  //vector<mpf_class> current_ukx  = compute_derivates(current_x,taylor_order);
  vector<mpf_class> current_ukx = start_ukx;
  mpf_class next_ux(0,precision);
  mpf_class next_dux(0,precision);
  mpf_class factorial(1,precision);
  mpf_class hk(1,precision);
  mpf_class h(next_x - current_x,precision);
  for(int i = 0; i < taylor_order.get_ui(); ++i)
  {
    next_ux = next_ux + (current_ukx.at(i)/factorial)*hk;
    next_dux = next_dux + (current_ukx.at(i+1)/factorial)*hk;
    factorial = factorial*mpf_class(i+1,precision);
    hk = hk*h;
  }
  next_ux = next_ux + (current_ukx.back()/factorial)*hk;
  vector<mpf_class> next_value = {next_ux,next_dux};
  return next_value;
};
#else
//do nothing；
#endif
