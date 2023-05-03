#ifndef __ORTHOGONAL_POLYNOMIALS__
#define __ORTHOGONAL_POLYNOMIALS__
#include "polynomial.h"
#include <functional>
//Family of Orthogonal Polynomials
class OrthogonalPolynomails
{
  protected:
    mpf_class interval[2];
    vector<polynomial> polynomials;
    vector<polynomial> derivate_polynomials;
    mpz_class highest_degree;
  public:
    OrthogonalPolynomails()
    {  
      highest_degree = mpz_class(1);
    };
    OrthogonalPolynomails(mpz_class _highestdegree)
    { 
      highest_degree = max(mpz_class(1),_highestdegree);
    };
    virtual ~OrthogonalPolynomails() = default;
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree)=0; /*< a_n,b_n,c_,n,d_n */
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree)=0;/*< p,q,r */
    virtual mpf_class weightfunction(mpf_class x) = 0;
    virtual polynomial getPm(mpz_class m)/*<get Pm polynomial of degree m*/
    {
      return polynomials.at(m.get_ui());
    }
    virtual polynomial getDerivatePm(mpz_class m)/*<get P'm polynomial*/
    {
      return derivate_polynomials.at(m.get_ui());
    }
    virtual void set_highest_degree(mpz_class degree)
    {
      if (highest_degree < degree) {
        polynomials.resize(degree.get_ui() + 1);
        derivate_polynomials.resize(degree.get_ui() + 1);
        for(int i = highest_degree.get_ui()+1; i <= degree.get_ui(); ++i)
        {
          vector<mpf_class> recur = coefficientsOfRecurrence(mpz_class(i-1));
          mpz_class degree_tmp(1);
          vector<mpf_class> coeff_tmp = {recur[1]/recur[0],recur[2]/recur[0]};
          polynomial tmp(coeff_tmp,degree_tmp);
          polynomial pi = tmp*polynomials.at(i-1) - (recur[3]/recur[0])*polynomials.at(i-2);
          polynomials.at(i) = pi;
          polynomial pi_d = tmp*derivate_polynomials.at(i-1) - (recur[3]/recur[0])*derivate_polynomials.at(i-2) + (recur[2]/recur[0])*polynomials.at(i-1);
          derivate_polynomials.at(i) = pi_d;
        }
        highest_degree = degree;
      }
      else {
        cout << "OrthogonalPolynomials already have degree "<< degree << endl;
        cout << "its highest_degree is " << highest_degree << endl;
      }
    };
  protected:
    virtual void initialize() = 0;
};


class LegendrePolys: public OrthogonalPolynomails
{
  public:
    LegendrePolys():OrthogonalPolynomails()
  {
    interval[0] = mpf_class(-1,precision);
    interval[1] = mpf_class(1,precision);
    initialize();
  };
    LegendrePolys(mpz_class _highestdegree)
      :OrthogonalPolynomails(_highestdegree)
    {
      interval[0] = mpf_class(-1,precision);
      interval[1] = mpf_class(1,precision);
      initialize(); 
    };
    ~LegendrePolys() = default;
    // degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree)
    {
      mpf_class a_n(degree+1,precision);
      mpf_class b_n(0,precision);
      mpf_class c_n(2*degree+1,precision);
      mpf_class d_n(degree,precision);
      vector<mpf_class> coefficients = {a_n,b_n,c_n,d_n};
      return coefficients;
    };
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree)
    {
      vector<mpf_class> p_coeff = {mpf_class(1,precision),mpf_class(0,precision),mpf_class(-1,precision)};
      mpz_class p_degree{2};
      polynomial p(p_coeff,p_degree);
      vector<mpf_class> q_coeff = {mpf_class(0,precision),mpf_class(-2,precision)};
      mpz_class q_degree{1};
      polynomial q(q_coeff,q_degree);
      vector<mpf_class> r_coeff = {mpf_class(degree*(degree+1),precision)};
      mpz_class r_degree{0};
      polynomial r(r_coeff,r_degree);
      vector<polynomial> coefficients = {p,q,r};
      return coefficients;
    };
    virtual mpf_class weightfunction(mpf_class x)
    {
      return mpf_class{1,precision};
    };
  private:
    //compute polynomial and derivate under highest_degree
    void initialize()
    {
      int psize = highest_degree.get_ui() + 1;
      polynomials.resize(psize);
      derivate_polynomials.resize(psize);
      //p0 = 1
      mpz_class degree_p0(0);
      vector<mpf_class> coeff_p0 = {mpf_class(1,precision)};
      polynomial p0(coeff_p0,degree_p0);
      polynomials.at(0) = p0;
      polynomial p0_d = derivate(p0);
      derivate_polynomials.at(0) = p0_d;
      //p1 = x
      mpz_class degree_p1(1);
      vector<mpf_class> coeff_p1 = {mpf_class(0,precision),mpf_class(1,precision)};
      polynomial p1(coeff_p1,degree_p1);
      polynomials.at(1) = p1;
      polynomial p1_d = derivate(p1);
      derivate_polynomials.at(1) = p1_d;
      for(int i = 2; i <= highest_degree; ++i)
      {
        vector<mpf_class> recur = coefficientsOfRecurrence(mpz_class(i-1));
        mpz_class degree_tmp(1);
        vector<mpf_class> coeff_tmp = {recur[1]/recur[0],recur[2]/recur[0]};
        polynomial tmp(coeff_tmp,degree_tmp);
        polynomial pi = tmp*polynomials.at(i-1) - (recur[3]/recur[0])*polynomials.at(i-2);
        polynomials.at(i) = pi;
        polynomial pi_d = tmp*derivate_polynomials.at(i-1) - (recur[3]/recur[0])*derivate_polynomials.at(i-2) + (recur[2]/recur[0])*polynomials.at(i-1);
        derivate_polynomials.at(i) = pi_d;
      }
    };
};

class HermitePolys: public OrthogonalPolynomails
{
  public:
    HermitePolys():OrthogonalPolynomails()
  {
    interval[0] = mpf_class(-infinity);
    interval[1] = mpf_class(infinity);
    initialize();
  };
    HermitePolys(mpz_class _highestdegree)
      :OrthogonalPolynomails(_highestdegree)
    {
      interval[0] = mpf_class(-infinity);
      interval[1] = mpf_class(infinity);
      initialize(); 
    };
    ~HermitePolys() = default;
    // degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree)
    {
      mpf_class a_n(1,precision);
      mpf_class b_n(0,precision);
      mpf_class c_n(2,precision);
      mpf_class d_n(2*degree,precision);
      vector<mpf_class> coefficients = {a_n,b_n,c_n,d_n};
      return coefficients;
    };
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree)
    {
      vector<mpf_class> p_coeff = {mpf_class(1,precision)};
      mpz_class p_degree{0};
      polynomial p(p_coeff,p_degree);
      vector<mpf_class> q_coeff = {mpf_class(0,precision),mpf_class(-2,precision)};
      mpz_class q_degree{1};
      polynomial q(q_coeff,q_degree);
      vector<mpf_class> r_coeff = {mpf_class(2*degree,precision)};
      mpz_class r_degree{0};
      polynomial r(r_coeff,r_degree);
      vector<polynomial> coefficients = {p,q,r};
      return coefficients;
    };
    virtual mpf_class weightfunction(mpf_class x)
    {
      mpf_class exception(-mpf_class_pow_ui(x, 2));
      double res = exp(exception.get_d());
      return mpf_class(res,precision);
    };
  private:
    //compute polynomial and derivate under highest_degree
    void initialize()
    {
      int psize = highest_degree.get_ui() + 1;
      polynomials.resize(psize);
      derivate_polynomials.resize(psize);
      //p0 = 1
      mpz_class degree_p0(0);
      vector<mpf_class> coeff_p0 = {mpf_class(1,precision)};
      polynomial p0(coeff_p0,degree_p0);
      polynomials.at(0) = p0;
      polynomial p0_d = derivate(p0);
      derivate_polynomials.at(0) = p0_d;
      //p1 = x
      mpz_class degree_p1(1);
      vector<mpf_class> coeff_p1 = {mpf_class(0,precision),mpf_class(2,256)};
      polynomial p1(coeff_p1,degree_p1);
      polynomials.at(1) = p1;
      polynomial p1_d = derivate(p1);
      derivate_polynomials.at(1) = p1_d;
      for(int i = 2; i <= highest_degree; ++i)
      {
        vector<mpf_class> recur = coefficientsOfRecurrence(mpz_class(i-1));
        mpz_class degree_tmp(1);
        vector<mpf_class> coeff_tmp = {recur[1]/recur[0],recur[2]/recur[0]};
        polynomial tmp(coeff_tmp,degree_tmp);
        polynomial pi = tmp*polynomials.at(i-1) - (recur[3]/recur[0])*polynomials.at(i-2);
        polynomials.at(i) = pi;
        polynomial pi_d = tmp*derivate_polynomials.at(i-1) - (recur[3]/recur[0])*derivate_polynomials.at(i-2) + (recur[2]/recur[0])*polynomials.at(i-1);
        derivate_polynomials.at(i) = pi_d;
      }
    };

};

class LaguerrePolys: public OrthogonalPolynomails
{
  public:
    LaguerrePolys():OrthogonalPolynomails()
  {
    interval[0] = mpf_class(0,precision);
    interval[1] = mpf_class(infinity,precision);
    initialize();
  };
    LaguerrePolys(mpz_class _highestdegree)
      :OrthogonalPolynomails(_highestdegree)
    {
      interval[0] = mpf_class(0,precision);
      interval[1] = mpf_class(infinity,precision);
      initialize(); 
    };
    ~LaguerrePolys() = default;
    // degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
    virtual vector<mpf_class> coefficientsOfRecurrence(mpz_class degree)
    {
      mpf_class a_n(degree+1,precision);
      mpf_class b_n(2*degree+1,precision);
      mpf_class c_n(-1,precision);
      mpf_class d_n(degree,precision);
      vector<mpf_class> coefficients = {a_n,b_n,c_n,d_n};
      return coefficients;
    };
    virtual vector<polynomial> coefficientsOfODE(mpz_class degree)
    {
      vector<mpf_class> p_coeff = {mpf_class(0),mpf_class(1)};
      mpz_class p_degree{1};
      polynomial p(p_coeff,p_degree);
      vector<mpf_class> q_coeff = {mpf_class(1,precision),mpf_class(-1,precision)};
      mpz_class q_degree{1};
      polynomial q(q_coeff,q_degree);
      vector<mpf_class> r_coeff = {mpf_class(degree,precision)};
      mpz_class r_degree{0};
      polynomial r(r_coeff,r_degree);
      vector<polynomial> coefficients = {p,q,r};
      return coefficients;
    };
    virtual mpf_class weightfunction(mpf_class x)
    {
      double res = exp(-x.get_d());
      return mpf_class(res,precision);
    };
  private:
    //compute polynomial and derivate under highest_degree
    void initialize()
    {
      int psize = highest_degree.get_ui() + 1;
      polynomials.resize(psize);
      derivate_polynomials.resize(psize);
      //p0 = 1
      mpz_class degree_p0(0);
      vector<mpf_class> coeff_p0 = {mpf_class(1,precision)};
      polynomial p0(coeff_p0,degree_p0);
      polynomials.at(0) = p0;
      polynomial p0_d = derivate(p0);
      derivate_polynomials.at(0) = p0_d;
      //p1 = x
      mpz_class degree_p1(1);
      vector<mpf_class> coeff_p1 = {mpf_class(1,precision),mpf_class(-1,precision)};
      polynomial p1(coeff_p1,degree_p1);
      polynomials.at(1) = p1;
      polynomial p1_d = derivate(p1);
      derivate_polynomials.at(1) = p1_d;
      for(int i = 2; i <= highest_degree; ++i)
      {
        vector<mpf_class> recur = coefficientsOfRecurrence(mpz_class(i-1));
        mpz_class degree_tmp(1);
        vector<mpf_class> coeff_tmp = {recur[1]/recur[0],recur[2]/recur[0]};
        polynomial tmp(coeff_tmp,degree_tmp);
        polynomial pi = tmp*polynomials.at(i-1) - (recur[3]/recur[0])*polynomials.at(i-2);
        polynomials.at(i) = pi;
        polynomial pi_d = tmp*derivate_polynomials.at(i-1) - (recur[3]/recur[0])*derivate_polynomials.at(i-2) + (recur[2]/recur[0])*polynomials.at(i-1);
        derivate_polynomials.at(i) = pi_d;
      }
    };

};
#else
//do nothing
#endif
