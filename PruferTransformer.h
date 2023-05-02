#ifndef __PRUFERTRANSFORMER__
#define __PRUFERTRANSFORMER__

#include "polynomial.h"
#include <climits>
#include <cstdint>
#include <functional>
#include <vector>
using namespace std;

const mpf_class pi(4.0*atan(1.0),256);
const mpf_class eps = numeric_limits<double>::epsilon();

class PruferTransformer
{
  private:
    function<mpf_class(mpf_class)> theta_x;
    function<mpf_class(mpf_class,mpf_class)>  dx_dtheta;
  public:
    PruferTransformer() = default;
    PruferTransformer(polynomial& u, vector<polynomial>& coeff_ODE);
    ~PruferTransformer() = default;
    function<mpf_class(mpf_class)> get_theta_x();
    function<mpf_class(mpf_class,mpf_class)> get_dx_dtheta();
    mpf_class theta_x_value(mpf_class x);
    mpf_class dx_dtheta_value(mpf_class theta,mpf_class x);
};

PruferTransformer::PruferTransformer(polynomial& u, vector<polynomial>& coeff_ODE)
{
  polynomial derivate_u = derivate(u);
  polynomial p(coeff_ODE.at(0));
  polynomial q(coeff_ODE.at(1));
  polynomial r(coeff_ODE.at(2));
  theta_x = [=](mpf_class x) mutable{
    mpf_class res;
    if (abs(r(x)) < eps || abs(u(x)) < eps) {
      res = pi/2.0;
    }
    else {
      mpf_class part1 = sqrt(p(x)/r(x));
      mpf_class part2 = derivate_u(x)/u(x);
      mpf_class mul = part1*part2;
      res = atan(mul.get_d());
    }

    return res;
  };
  dx_dtheta = [=](mpf_class theta,mpf_class x) mutable{
    polynomial derivate_p = derivate(p);
    polynomial derivate_r = derivate(r);
    polynomial tmp1 = derivate_r*p -  derivate_p*r + 2*r*q;
    mpf_class rx = r(x), px = p(x),tmp1x = tmp1(x);
    mpf_class part1,part2,part3,part4;
    mpf_class denominator(0);
    if (abs(rx) < eps) {
      part1 = 0;
    }
    else if (abs(px) < eps) {
      part1 = INT64_MAX;
      denominator = INT64_MAX;
    }
    else {
      part1 = sqrt(r(x)/p(x));
    }
    part2 = tmp1x*sin(2*theta.get_d());
    part3 = 4*rx*px;
    if (abs(part2) < eps) {
      part4 = 0;
    }
    else if (abs(part3) < eps) {
      part4 = INT64_MAX;
      denominator = INT64_MAX;
    }
    else {
      part4 = part2/part3;
    }

    if(denominator == 0)
    {
      denominator = part1 + part4;
    }
    mpf_class res;
    if (abs(denominator) < eps) {
      res = -INT64_MAX;
    }
    else if (denominator == INT64_MAX) {
      res = 0;
    }
    else {
      res = -mpf_class(1.0)/denominator;
    }
    return res;
  };
}

function<mpf_class(mpf_class)> PruferTransformer::get_theta_x()
{
  return theta_x;
}

function<mpf_class(mpf_class,mpf_class)> PruferTransformer::get_dx_dtheta()
{
  return dx_dtheta;
}

mpf_class PruferTransformer::theta_x_value(mpf_class x)
{
  return theta_x(x);
};

mpf_class PruferTransformer::dx_dtheta_value(mpf_class theta,mpf_class x)
{ 
  return dx_dtheta(theta,x);
};
#else
//do nothing
#endif
