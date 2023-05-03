#include <vector>
#include "GaussianPoint.h"
#include "gmptools.h"
#include <algorithm>
#include <iostream>
#include "GaussIntergalTableGenerator.h"
using namespace std;
void test_copy()
{
  vector<mpf_class> half_odd = {4,5,6,7,8};
  vector<mpf_class> odd(9);
  reverse_copy(half_odd.begin(), half_odd.end(), odd.begin());
  copy(half_odd.begin()+1, half_odd.end(), odd.begin()+4+1);
  for (mpf_class m : odd)
  {
    cout << m << " ";
  }
  cout << endl;
  vector<mpf_class> half_even = {5,6,7,8,9};
  vector<mpf_class> even(10);
  reverse_copy(half_even.begin(), half_even.end(), even.begin());
  copy(half_even.begin(), half_even.end(), even.begin()+4+1);
  for (mpf_class m : even) {
    cout << m << " ";
  }
  cout << endl;
}

void test_Gauss_Legendre_3()
{
  mpz_class highest_prec = 5;
  OrthogonalPolynomails* Pn = new LegendrePolys(3);
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (auto gauss : gausstable)
  {
    cout << gauss;
  }
  //delete Pn;
}

void test_Gauss_Legendre_4()
{
  mpz_class highest_prec = 7;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LegendrePolys(2);
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (auto gauss : gausstable)
  {
    cout << gauss;
  }
  //delete Pn;
}
int main()
{
  test_copy();
  test_Gauss_Legendre_3();
  test_Gauss_Legendre_4();
  return 0;
}
