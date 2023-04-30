#ifndef __ORTHOGONAL_POLYNOMIALS__
#define __ORTHOGONAL_POLYNOMIALS__
#include "polynomial.h"
#include <functional>
#include <gmpxx.h>
//Family of Orthogonal Polynomials
class OrthogonalPolynomails
{
  private:
    mpf_class interval[2];
    //function<mpf_class(mpf_class)> weightfunction;
    vector<polynomial> polynomials;
    mpz_class highestdegree;
  public:
    OrthogonalPolynomails()
    {  
      interval[0] = mpf_class(-1.0);
      interval[1] = mpf_class(1.0);
      highestdegree = mpz_class(0);
    };
    OrthogonalPolynomails(mpf_class _interval[],mpz_class _highestdegree)
    {
      for(int i = 0; i < 2; ++i) {
        interval[i] = _interval[i];
      }
      highestdegree = _highestdegree;
    };
    virtual ~OrthogonalPolynomails();
    virtual vector<mpz_class> coefficientsOfRecurrence(mpz_class degree)=0; /*< a_n,b_n,c_,n,d_n */
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree)=0;/*< p,q,r */
    virtual mpf_class weightfunction(mpf_class x) = 0;
};


class LegendrePolys: public OrthogonalPolynomails
{
  public:
    LegendrePolys():OrthogonalPolynomails(){};
    LegendrePolys(mpf_class _interval[],mpz_class _highestdegree)
      :OrthogonalPolynomails(_interval,_highestdegree)
    {};
    ~LegendrePolys();
    // degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
    virtual vector<mpz_class> coefficientsOfRecurrence(mpz_class degree)
    {
      mpz_class a_n(degree+1);
      mpz_class b_n(0);
      mpz_class c_n(2*degree+1);
      mpz_class d_n(degree);
      vector<mpz_class> coefficients = {a_n,b_n,c_n,d_n};
      return coefficients;
    };
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree)
    {
      vector<mpf_class> p_coeff = {1,0,-1};
      mpz_class p_degree{2};
      polynomial p(p_coeff,p_degree);
      vector<mpf_class> q_coeff = {0,1};
      mpz_class q_degree{1};
      polynomial q(q_coeff,q_degree);
      vector<mpf_class> r_coeff = {degree*(degree+1)};
      mpz_class r_degree{0};
      polynomial r(r_coeff,r_degree);
      vector<polynomial> coefficients = {p,q,r};
      return coefficients;
    };
    virtual mpf_class weightfunction(mpf_class x)
    {
      return mpf_class{1};
    };
};

class HermitePolys: public OrthogonalPolynomails
{

};

class LaguerrePolys: public OrthogonalPolynomails
{

};
#else
//do nothing
#endif
