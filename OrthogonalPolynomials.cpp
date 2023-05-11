#include "OrthogonalPolynomials.h"
#include "iostream"
#include "math.h"
using namespace std;

//正交多项式族类
OrthogonalPolynomials::OrthogonalPolynomials()
{  
  highest_degree = mpz_class(1);
};
OrthogonalPolynomials::OrthogonalPolynomials(mpz_class _highestdegree)
{ 
  highest_degree = max(mpz_class(1),_highestdegree);
};
/*get Pm polynomial of degree m*/
polynomial OrthogonalPolynomials::getPm(mpz_class m)
{
  return polynomials.at(m.get_ui());
}
/*get derivate polynomial of Pm */
polynomial OrthogonalPolynomials::getDerivatePm(mpz_class m)
{
  return derivate_polynomials.at(m.get_ui());
}
void OrthogonalPolynomials::set_highest_degree(mpz_class degree)
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
vector<mpf_class> OrthogonalPolynomials::getQuadratureWeight(mpz_class degree, vector<mpf_class> roots)
{
  int len = roots.size();
  vector<mpf_class> weights(len,mpf_class(0,2*precision));
  for (int i = 0; i < len; ++i) {
    mpf_class root(0,2*precision);
    root = roots.at(i);
    weights.at(i) = getQuadratureWeight(degree,root);
  }
  return weights;
};


//勒让德多项式


LegendrePolys::LegendrePolys():OrthogonalPolynomials()
{
  interval[0] = mpf_class(-1,precision);
  interval[1] = mpf_class(1,precision);
  initialize();
};
LegendrePolys::LegendrePolys(mpz_class _highestdegree)
  :OrthogonalPolynomials(_highestdegree)
{
  interval[0] = mpf_class(-1,precision);
  interval[1] = mpf_class(1,precision);
  initialize(); 
};
// degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
vector<mpf_class> LegendrePolys::coefficientsOfRecurrence(mpz_class degree)
{
  mpf_class a_n(degree+1,precision);
  mpf_class b_n(0,precision);
  mpf_class c_n(2*degree+1,precision);
  mpf_class d_n(degree,precision);
  vector<mpf_class> coefficients = {a_n,b_n,c_n,d_n};
  return coefficients;
};
vector<polynomial> LegendrePolys::coefficientsOfODE(mpz_class degree)
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
mpf_class LegendrePolys::weightfunction(mpf_class x)
{
  return mpf_class{1,precision};
};
mpf_class LegendrePolys::getQuadratureWeight(mpz_class degree, mpf_class root)
{
  polynomial dp = derivate_polynomials.at(degree.get_ui());
  mpf_class weight(0,precision);
  mpf_class dpx = dp(root);

  mpf_class denominator(0,precision);
  denominator = (1.0-mpf_class_pow_ui(root, 2))*mpf_class_pow_ui(dpx, 2);
  weight = mpf_class(2.0,precision)/denominator;
  return weight;
};

//compute polynomial and derivate under highest_degree
void LegendrePolys::initialize()
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


// 赫米特多项式
HermitePolys::HermitePolys():OrthogonalPolynomials()
{
  interval[0] = mpf_class(-infinity);
  interval[1] = mpf_class(infinity);
  initialize();
};
HermitePolys::HermitePolys(mpz_class _highestdegree)
  :OrthogonalPolynomials(_highestdegree)
{
  interval[0] = mpf_class(-infinity);
  interval[1] = mpf_class(infinity);
  initialize(); 
};
// degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
vector<mpf_class> HermitePolys::coefficientsOfRecurrence(mpz_class degree)
{
  mpf_class a_n(1,precision);
  mpf_class b_n(0,precision);
  mpf_class c_n(2,precision);
  mpf_class d_n(2*degree,precision);
  vector<mpf_class> coefficients = {a_n,b_n,c_n,d_n};
  return coefficients;
};
vector<polynomial> HermitePolys::coefficientsOfODE(mpz_class degree)
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
mpf_class HermitePolys::weightfunction(mpf_class x)
{
  mpf_class exception(-mpf_class_pow_ui(x, 2));
  double res = exp(exception.get_d());
  return mpf_class(res,precision);
};
mpf_class HermitePolys::getQuadratureWeight(mpz_class degree, mpf_class root)
{
  int order = degree.get_ui();
  mpf_class weight(0,precision);
  //TODO: write code
  mpf_class sqrt_pi(0,precision);
  sqrt_pi = sqrt(pi);
  mpf_class fact = mpf_class_pow_ui(mpf_class(2.0,precision), order+1)*factorial(degree);
  polynomial dh = derivate_polynomials.at(order);
  mpf_class dhx = dh(root);
  weight = (sqrt_pi*fact)/(dhx*dhx);
  return weight;
}
//compute polynomial and derivate under highest_degree
void HermitePolys::initialize()
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


//拉盖尔多项式

LaguerrePolys::LaguerrePolys():OrthogonalPolynomials()
{
  interval[0] = mpf_class(0,precision);
  interval[1] = mpf_class(infinity,precision);
  initialize();
};
LaguerrePolys::LaguerrePolys(mpz_class _highestdegree)
  :OrthogonalPolynomials(_highestdegree)
{
  interval[0] = mpf_class(0,precision);
  interval[1] = mpf_class(infinity,precision);
  initialize(); 
};
// degree == n,a_{n}*P_{n+1} = (b_{n}+c_{n}*x)*P_{n}-d_{n}P_{n-1}
vector<mpf_class> LaguerrePolys::coefficientsOfRecurrence(mpz_class degree)
{
  mpf_class a_n(degree+1,precision);
  mpf_class b_n(2*degree+1,precision);
  mpf_class c_n(-1,precision);
  mpf_class d_n(degree,precision);
  vector<mpf_class> coefficients = {a_n,b_n,c_n,d_n};
  return coefficients;
};
vector<polynomial> LaguerrePolys::coefficientsOfODE(mpz_class degree)
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
mpf_class LaguerrePolys::weightfunction(mpf_class x)
{
  double res = exp(-x.get_d());
  return mpf_class(res,precision);
};
mpf_class LaguerrePolys::getQuadratureWeight(mpz_class degree, mpf_class root)
{
  mpf_class weight(0,precision);
  int order = degree.get_ui();
  polynomial dl = derivate_polynomials.at(order);
  mpf_class dlx1(0,2*precision);
  dlx1 = dl(root);
  mpf_class dlx2(0,2*precision);
  dlx2 = dlx1*dlx1;
  mpf_class root1(0,2*precision);
  root1 = root;
  mpf_class denominator(0,2*precision);
  denominator = root1*dlx2;
  weight = mpf_class(1.0,2*precision)/denominator;
  return weight;
  // another way to compute weight
  // mpf_class weight(0,precision);
  // int order = degree.get_ui();
  // mpf_class dlx(0,8*precision);
  // dlx = polynomials.at(order-1)(root);
  // mpf_class dlx2(0,4*precision);
  // dlx2 = dlx*dlx;
  // mpf_class n2((degree)*(degree),4*precision);
  // mpf_class denominator(0,2*precision);
  // denominator = n2*dlx2;
  // weight = mpf_class(root,2*precision)/denominator;
  // return weight;
};
//compute polynomial and derivate under highest_degree
void LaguerrePolys::initialize()
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
    // polynomial pi_d = derivate(pi);
    derivate_polynomials.at(i) = pi_d;
  }
};
