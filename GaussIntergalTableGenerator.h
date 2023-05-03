#ifndef __GAUSSIANINTERGRALTABLEGENERATOR__
#define __GAUSSIANINTERGRALTABLEGENERATOR__

#include "GaussianPoint.h"
#include "OrthogonalPolynomials.h"
#include "PruferTransformer.h"
#include "RungeKutta.h"
#include "NewtonMethod.h"
#include <algorithm>
#include <gmp.h>
#include <memory>
class GaussIntegralTableGenerator 
{
  private:
    mpz_class highest_prec;
    unique_ptr<OrthogonalPolynomails> ortho_polys;
    vector<GaussianPoint1D> table;
  public:
    GaussIntegralTableGenerator();
    GaussIntegralTableGenerator(mpz_class _hightest_prec, OrthogonalPolynomails* _ortho_polys);
    const vector<GaussianPoint1D>& compute_gaussian_table();
    GaussianPoint1D compute_gaussian_table(mpz_class spec_prec);
  private:
    mpf_class compute_first_root(mpz_class poly_order);
    mpf_class compute_subsequent_root(mpz_class poly_order,mpf_class current_x);
};

GaussIntegralTableGenerator::GaussIntegralTableGenerator()
  :highest_prec(1)
{
  ortho_polys = make_unique<LegendrePolys>(1);
};

GaussIntegralTableGenerator::GaussIntegralTableGenerator(mpz_class _highest_prec, OrthogonalPolynomails* _ortho_polys)
{
  highest_prec = max(mpz_class(1),_highest_prec);
  unique_ptr<OrthogonalPolynomails> ptr(_ortho_polys);
  ortho_polys = std::move(ptr);
}

const vector<GaussianPoint1D>& GaussIntegralTableGenerator::compute_gaussian_table()
{
  mpz_class poly_order = ceil((highest_prec.get_d()+1)/2.0);
  table.resize(poly_order.get_ui() + 1);
  ortho_polys->set_highest_degree(poly_order);
  for (mpz_class order = 1; order <= poly_order; ++order) {
    GaussianPoint1D gauss(compute_gaussian_table(order));
    table.at(order.get_ui()) = gauss;
  }
  return table;
};

GaussianPoint1D GaussIntegralTableGenerator::compute_gaussian_table(mpz_class spec_order)
{
  mpz_class point_num = spec_order;
  mpz_class is_odd = point_num%2;
  ortho_polys->set_highest_degree(spec_order);
  // 如果是Legendre多项式
  polynomial p = ortho_polys->getPm(spec_order);
  polynomial dp = ortho_polys->getDerivatePm(spec_order);
  mpz_class x_num = ceil(point_num.get_d()/2.0) - 1;
  // 奇偶性，先求大于等于0的一侧，然后对称。
  vector<mpf_class> half_gauss_point{x_num+1,mpf_class(0,precision)};
  mpf_class x0 = compute_first_root(spec_order);
  half_gauss_point.at(0) = x0;
  //if n is odd, (n-1)/2 = (n+1)/2 - 1;
  //if n is even ,(n)/2;
  for(int i = 1; i <= x_num.get_ui(); ++i) {
    mpf_class next_x = compute_subsequent_root(spec_order,half_gauss_point.at(i-1));
    half_gauss_point.at(i) = next_x;
  }
  //TODO:计算积分权重
  vector<mpf_class> half_weight{x_num+1,mpf_class(0,precision)};
  vector<mpf_class> weight{point_num,mpf_class(0,precision)};

  vector<mpf_class> gauss_point{point_num,mpf_class(0,precision)};
  if (is_odd !=0) {
    reverse_copy(half_gauss_point.begin(), half_gauss_point.end(),gauss_point.begin());
    copy(half_gauss_point.begin()+1,half_gauss_point.end(),gauss_point.begin() + x_num.get_ui()+1);
  } 
  else {
    reverse_copy(half_gauss_point.begin(), half_gauss_point.end(),gauss_point.begin());
    copy(half_gauss_point.begin(),half_gauss_point.end(),gauss_point.begin() + x_num.get_ui()+1);
 
  }

  GaussianPoint1D gauss(point_num,gauss_point,weight);
  return gauss;
};

mpf_class GaussIntegralTableGenerator::compute_first_root(mpz_class poly_order)
{
  //如果是Legendre多项式
  mpz_class is_odd = poly_order%2;
  // 对于奇数阶,x=0即为第一个根。
  if (is_odd != 0) {
    return mpf_class{0,precision};
  }
  //对于偶数阶，x=0处是个极值，使用RK和Newton法求第一个根
  else {
    polynomial p = ortho_polys->getPm(poly_order);
    vector<polynomial> coeff_ODE = ortho_polys->coefficientsOfODE(poly_order);
    PruferTransformer prufer(p,coeff_ODE);
    auto dx_dtheta = prufer.get_dx_dtheta();
    vector<mpf_class> initial_value = {mpf_class(0,precision),prufer.theta_x_value(0)};
    // use RungeKutta4 to approximate first root x0
    RungeKutta4 RK4(initial_value,dx_dtheta);
    mpf_class h = -1e-3;
    mpf_class end_x = -pi/2.0;
    mpf_class end_y = RK4.compute(h, end_x);
    polynomial dp = ortho_polys->getDerivatePm(poly_order);
    vector<polynomial> initial_function = {p,dp};
    // improve the precision of x0 via Newton's method
    NewtonMethod newton(end_y,initial_function,coeff_ODE);
    mpf_class first_root = newton.compute_root_with_taylor(120,80); 
    return first_root;
  }
};

mpf_class GaussIntegralTableGenerator::compute_subsequent_root(mpz_class poly_order,mpf_class current_x)
{
    polynomial p = ortho_polys->getPm(poly_order);
    vector<polynomial> coeff_ODE = ortho_polys->coefficientsOfODE(poly_order);
    PruferTransformer prufer(p,coeff_ODE);
    auto dx_dtheta = prufer.get_dx_dtheta();
    vector<mpf_class> initial_value = {mpf_class(0.5*pi,precision),current_x};
    // use RungeKutta4 to approximate first root x0
    RungeKutta4 RK4(initial_value,dx_dtheta);
    mpf_class h = -1e-3;
    mpf_class end_x = -pi/2.0;
    mpf_class end_y = RK4.compute(h, end_x);
    vector<polynomial> initial_function = {p,ortho_polys->getDerivatePm(poly_order)};
    // improve the precision of x0 via Newton's method
    NewtonMethod newton(end_y,initial_function,coeff_ODE);
    mpf_class next_root = newton.compute_root_with_taylor(120,80);
    return next_root;
};
#else
//do nothing
#endif
