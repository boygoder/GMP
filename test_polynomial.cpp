#include "polynomial.h"
#include <gmpxx.h>
#include <iostream>
#include <assert.h>
#include <iomanip>
void test_default_construct()
{
  cout <<"test default constructor of polynomial:\n";
  polynomial p1;
  assert(p1.getDegree()==0);
  cout << "p1: " << p1 << endl;
}

void test_construct()
{
  cout << "test construct of polynomial:\n";
  mpz_class degree(2);
  vector<mpf_class> coeff = {1.0,2.0,1.0};
  polynomial p1(coeff,degree);
  assert(p1.getDegree()==degree);
  cout << "p1: " <<p1 << endl;
}

void test_add_operator()
{
  cout << "test add operator of polynomial:\n";
  mpz_class degree_p1(2);
  vector<mpf_class> coeff_p1 = {1,2,1};
  polynomial p1(coeff_p1,degree_p1);
  mpz_class degree_p2(1);
  vector<mpf_class> coeff_p2 = {-2,-1};
  polynomial p2(coeff_p2,degree_p2);
  polynomial p3 = p1+p2;
  cout << p1 << p2 << endl;
  cout << "p1+p2=" << endl;
  cout << p3; 
  assert(p3.getDegree() == 2);
  vector<mpf_class> coeff_p3 = {-1,1,1};
  assert(p3.getCoeff() == coeff_p3);
}

void test_sub_operator()
{
  cout << "test sub operator of polynomial\n";
  mpz_class degree_p1(2);
  vector<mpf_class> coeff_p1 = {1,1,-1};
  polynomial p1(coeff_p1,degree_p1);
  mpz_class degree_p2(2);
  vector<mpf_class> coeff_p2 = {1,2,1};
  polynomial p2(coeff_p2,degree_p2);
  polynomial p3 = p1 - p2;
  vector<mpf_class> coeff_p3 = {0,-1,-2};
  cout << p1 << p2;
  cout <<"p1-p2=\n";
  cout << p3 << endl;
  assert(p3.getDegree() == 2);
  assert(p3.getCoeff() == coeff_p3);

}

void test_multi_digit()
{
  cout <<"test multiply digit of polynomial\n";
  mpf_class digit(2.5);
  mpz_class degree_p1(2);
  vector<mpf_class> coeff_p1 = {1,2,1};
  polynomial p1(coeff_p1,degree_p1);
  polynomial p2 = digit*p1;
  vector<mpf_class> coeff_p2 = {2.5,5.0,2.5};
  cout << p1;
  cout << digit <<"*p1 =\n";
  cout << p2;
  assert(p2.getDegree() == 2);
  assert(p2.getCoeff() == coeff_p2);
}

void test_multi_poly()
{
  cout << "test multiply polynomial of polynomial :\n";
  mpz_class degree_p1(1);
  vector<mpf_class> coeff_p1 = {2,1};
  polynomial p1(coeff_p1,degree_p1);
  mpz_class degree_p2(1);
  vector<mpf_class> coeff_p2 = {1,1};
  polynomial p2(coeff_p2,degree_p2);
  polynomial p3 = p1*p2;
  vector<mpf_class> coeff_p3 = {2,3,1};
  cout << p1 << p2;
  cout << "p1*p2=\n";
  cout << p3;
  assert(p3.getDegree() == 2);
  assert(p3.getCoeff() == coeff_p3);
}

void test_point_value()
{
  cout << "test pointValue of polynomial:\n";
  mpz_class degree_p1(2);
  vector<mpf_class> coeff_p1 = {1,2,1};
  polynomial p1(coeff_p1,degree_p1);
  cout << p1;
  cout <<"p1(-1) = "<< pointValue(mpf_class(-1),p1) <<"\n";
  cout <<"p1(-2) = "<< pointValue(mpf_class(-2),p1) <<"\n";
  cout <<"p1(0) = "<< pointValue(mpf_class(0),p1) <<"\n";
  assert(pointValue(mpf_class(-1), p1) == p1(-1));
  assert(pointValue(mpf_class(-2), p1) == p1(-2));
  assert(pointValue(mpf_class(0), p1) == p1(0));
}

void test_derivate()
{
  cout << "test derivate of polynomial:\n";
  mpz_class degree_p1(2);
  vector<mpf_class> coeff_p1 = {3,2,1};
  polynomial p1(coeff_p1,degree_p1); 
  polynomial p1_d = derivate(p1);
  vector<mpf_class> coeff_p1d = {2,2};
  cout << p1;
  cout <<"derivate of p1 is:\n";
  cout << p1_d;
  assert(p1_d.getDegree() == 1);
  assert(p1_d.getCoeff() == coeff_p1d);
  
  mpz_class degree_p2(0);
  vector<mpf_class> coeff_p2 = {1};
  polynomial p2 = polynomial(coeff_p2,degree_p2);
  polynomial p2_d = derivate(p2);
  vector<mpf_class> coeff_p2d = {0};
  cout << p2;
  cout << "derivate of p2 is:\n";
  cout << p2_d;
  assert(p2_d.getDegree() == 0);
  assert(p2_d.getCoeff() == coeff_p2d);
}

void test_brackets_operator()
{
  mpz_class degree_p1(2);
  vector<mpf_class> coeff_p1 = {1,2,1};
  polynomial p1(coeff_p1,degree_p1);
  cout << "p1(x=-1) = " << p1(mpf_class(-1)) << endl;
  cout << "p1(x=0) = " << p1(mpf_class(0)) << endl;
  cout << "p1(x=1) = " << p1(mpf_class(1)) << endl;
  cout << "p1(x=2) = " << p1(mpf_class(2)) << endl;
  polynomial p2(p1);
  mpf_class v2 = sqrt(p2(1)*p2(2));
  assert(p1(-1) == 0);
  assert(p1(0) == 1);
  assert(p1(1) == 4);
  assert(p1(2) == 9);
}
void test_new_pointvalue()
{
  cout << "test new pointValue:\n";
  cout << fixed << setprecision(80);
  mpz_class degree_p1(4);
  vector<mpf_class> coeff_p1 = {-4,3,-3,0,2};
  polynomial p1(coeff_p1,degree_p1);
  cout << "p1(-2) = " << p1(-2) << "\n";
  assert(p1(-2) == 10);
  cout << "p1'(-2) = " << p1.derivateValue(-2) << "\n";
  assert(p1.derivateValue(-2) == -49);
  cout << endl;
}
int main()
{
  test_default_construct();
  test_construct();
  test_brackets_operator();
  test_add_operator();
  test_sub_operator();
  test_multi_digit();
  test_multi_poly();
  test_point_value();
  test_derivate();
  test_new_pointvalue();
  return 0;
}

