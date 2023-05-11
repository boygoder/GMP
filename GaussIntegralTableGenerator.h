#ifndef __GAUSSIANINTERGRALTABLEGENERATOR__
#define __GAUSSIANINTERGRALTABLEGENERATOR__

#include "GaussianPoint.h"
#include "OrthogonalPolynomials.h"
#include <memory>
class GaussIntegralTableGenerator 
{
  private:
    mpz_class highest_prec;
    unique_ptr<OrthogonalPolynomials> ortho_polys;
    vector<GaussianPoint1D> table;
    mpz_class taylor_items = 120;
    mpz_class newton_eps = 80;
  public:
    GaussIntegralTableGenerator();
    GaussIntegralTableGenerator(mpz_class _hightest_prec, OrthogonalPolynomials* _ortho_polys);
    const vector<GaussianPoint1D>& compute_gaussian_table();
    GaussianPoint1D compute_gaussian_table(mpz_class spec_prec);
    void set_taylor_items(mpz_class _taylor_items);
    void set_newton_eps(mpz_class _newton_eps);
  private:
    mpf_class compute_first_root_symmetric(mpz_class poly_order);
    mpf_class compute_first_root_unsymmetric(mpz_class poly_order);
    mpf_class compute_subsequent_root(mpz_class poly_order,mpf_class current_x);
    mpf_class compute_subsequent_root_unsymmetric(mpz_class poly_order,mpf_class current_x);
    GaussianPoint1D compute_gaussian_symmetric(mpz_class spec_order);
    GaussianPoint1D compute_gaussian_unsymmetric(mpz_class spec_order);
};


#else
//do nothing
#endif
