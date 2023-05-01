#ifndef __PRUFERTRANSFORMER__
#define __PRUFERTRANSFORMER__

#include "polynomial.h"
#include <functional>
#include <vector>
using namespace std;
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
};

PruferTransformer::PruferTransformer(polynomial& u, vector<polynomial>& coeff_ODE)
{
  polynomial derivate_u = derivate(u);
  polynomial p(coeff_ODE.at(0));
  polynomial q(coeff_ODE.at(1));
  polynomial r(coeff_ODE.at(2));
  theta_x = [=](mpf_class x) mutable{
    mpf_class part1 = 1.0/(sqrt(r(x)*p(x)));
    mpf_class part2 = (p(x)*derivate_u(x))/(u(x));
    mpf_class mul = part1*part2;
    mpf_class res = atan(mul.get_d());
    return res;
  };
  dx_dtheta = [=](mpf_class theta,mpf_class x) mutable{
    polynomial derivate_p = derivate(p);
    polynomial derivate_r = derivate(u);
    mpf_class part1 = sqrt(r(x)/p(x));
    mpf_class part2 = (derivate_r(x)*p(x) - derivate_p(x)*r(x) + 2*r(x)*q(x))/(2*r(x)*p(x));
    mpf_class part3 = sin(2*theta.get_d())/2.0;
    mpf_class res = -1.0/(part1 + part2*part3);
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


#else
//do nothing
#endif
