#ifndef __ORTHOGONAL_POLYNOMIALS__
#define __ORTHOGONAL_POLYNOMIALS__
#include "gmptools.h"
#include "polynomial.h"
#include <functional>
//Family of Orthogonal Polynomials
class OrthogonalPolynomials
{
  protected:
    mpf_class interval[2];
    vector<polynomial> polynomials;
    vector<polynomial> derivate_polynomials;
    mpz_class highest_degree;
  public:
    OrthogonalPolynomials();
    OrthogonalPolynomials(mpz_class _highestdegree);
    virtual ~OrthogonalPolynomials() = default;
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree)=0; /*< a_n,b_n,c_,n,d_n */
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree)=0;/*< p,q,r */
    virtual mpf_class weightfunction(mpf_class x) = 0;
    /*get Pm polynomial of degree m*/
    virtual polynomial getPm(mpz_class m);
    /*get derivate polynomial of Pm */
    virtual polynomial getDerivatePm(mpz_class m);
    virtual void set_highest_degree(mpz_class degree);
    virtual vector<mpf_class> getQuadratureWeight(mpz_class degree, vector<mpf_class> roots);
    virtual mpf_class getQuadratureWeight(mpz_class degree, mpf_class root) = 0;
  protected:
    virtual void initialize() = 0;
};


class LegendrePolys: public OrthogonalPolynomials
{
  public:
    LegendrePolys();
    LegendrePolys(mpz_class _highestdegree);
    ~LegendrePolys() = default;
    // degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree);
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree);
    virtual mpf_class weightfunction(mpf_class x);
    virtual mpf_class getQuadratureWeight(mpz_class degree, mpf_class root);
  private:
    //compute polynomial and derivate under highest_degree
    virtual void initialize();
};

class HermitePolys: public OrthogonalPolynomials
{
  public:
    HermitePolys();
    HermitePolys(mpz_class _highestdegree);
    ~HermitePolys() = default;
    // degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree);
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree);
    virtual mpf_class weightfunction(mpf_class x);
    virtual mpf_class getQuadratureWeight(mpz_class degree, mpf_class root);
  private:
    //compute polynomial and derivate under highest_degree
    virtual void initialize();

};

class LaguerrePolys: public OrthogonalPolynomials
{
  public:
    LaguerrePolys();
    LaguerrePolys(mpz_class _highestdegree);
    ~LaguerrePolys() = default;
    // degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree);
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree);
    virtual mpf_class weightfunction(mpf_class x);
    virtual mpf_class getQuadratureWeight(mpz_class degree, mpf_class root);
  private:
    //compute polynomial and derivate under highest_degree
    virtual void initialize();

};
#else
//do nothing
#endif
