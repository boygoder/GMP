#include "PruferTransformer.h"
#include "math.h"
PruferTransformer::PruferTransformer(polynomial& u, vector<polynomial>& coeff_ODE)
{
  polynomial derivate_u = derivate(u);
  polynomial p(coeff_ODE.at(0));
  polynomial q(coeff_ODE.at(1));
  polynomial r(coeff_ODE.at(2));
  theta_x = [=](mpf_class x) mutable{
    mpf_class res(0,precision);
    if (abs(r(x)) < zero_eps|| abs(u(x)) < zero_eps) {
      res = pi/mpf_class(2.0,precision);
    }
    else {
      mpf_class part1 = sqrt(p(x)/r(x));
      mpf_class part2 = derivate_u(x)/u(x);
      mpf_class mul = part1*part2;
      //we use atan function in std, which reduces the precision
      // we only use this function for RungeKutta.
      res = atan(mul.get_d());
    }

    return res;
  };
  dx_dtheta = [=](mpf_class theta,mpf_class x) mutable{
    polynomial derivate_p = derivate(p);
    polynomial derivate_r = derivate(r);
    polynomial tmp1 = derivate_r*p -  derivate_p*r + 2*r*q;
    mpf_class rx = r(x), px = p(x),tmp1x = tmp1(x);
    mpf_class part1(0,precision),part2(0,precision);
    mpf_class part3(0,precision),part4(0,precision);
    mpf_class denominator(0,precision);
    if (abs(rx) < zero_eps) {
      part1 = 0;
    }
    else if (abs(px) < zero_eps) {
      part1 = infinity;
      // denominator = infinity;
    }
    else {
      part1 = sqrt(r(x)/p(x));
    }
    part2 = tmp1x*sin(2*theta.get_d());
    part3 = 4*rx*px;
    if (abs(part2) < zero_eps) {
      part4 = 0;
    }
    else if (abs(part3) < zero_eps) {
      part4 = infinity;
      // denominator = infinity;
    }
    else {
      part4 = part2/part3;
    }
    // if(denominator == 0)
    // {
      denominator = part1 + part4;
    // }
    mpf_class res;
    if (abs(denominator) < zero_eps) {
      res = -infinity;
    }
    else if (denominator >= infinity) {
      res = 0;
    }
    else {
      res = -mpf_class(1,precision)/denominator;
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