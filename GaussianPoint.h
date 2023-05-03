#ifndef __GAUSSIANINTERGRALTABLEGENERATOR__
#define __GAUSSIANINTERGRALTABLEGENERATOR__


#include "OrthogonalPolynomials.h"
#include "RungeKutta.h"
#include "NewtonMethod.h"
#include <memory>
class GaussIntegralTableGenerator 
{
  private:
    mpz_class alg_prec;
    unique_ptr<OrthogonalPolynomails> ortho_polys;
  public:
    GaussIntegralTableGenerator();
    GaussIntegralTableGenerator(mpz_class _alg_prec, OrthogonalPolynomails* _ortho_polys);
    void compute_gaussian_table();
  private:
    mpf_class compute_first_root();
    mpf_class compute_subsequent_root();
};


GaussIntegralTableGenerator::GaussIntegralTableGenerator()
{
  alg_prec = 1;
  ortho_polys = make_unique<LegendrePolys>(2);
}
GaussIntegralTableGenerator::GaussIntegralTableGenerator(mpz_class _alg_prec, OrthogonalPolynomails* _ortho_polys)
{
  alg_prec = max(mpz_class(1),_alg_prec);
  unique_ptr<OrthogonalPolynomails> ptr(_ortho_polys);
  ortho_polys = std::move(ptr);
}
#else
//do nothing
#endif
