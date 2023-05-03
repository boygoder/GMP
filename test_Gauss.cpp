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

void test_transform()
{
  vector<mpf_class> half_odd = {4,5,6,7,8};
  vector<mpf_class> odd(9,mpf_class(0,precision));
  reverse_copy(half_odd.begin()+1, half_odd.end(), odd.begin());
  copy(half_odd.begin(), half_odd.end(), odd.begin()+4);
  std::transform(odd.begin(), odd.begin() + 4, odd.begin(), 
        [](mpf_class i) { 
          return mpf_class(-1.0*i,precision); 
        });
  for (mpf_class m : odd)
  {
    cout << m << "\n";
  }
  cout << endl;
  vector<mpf_class> half_even = {5,6,7,8,9};
  vector<mpf_class> even(10,mpf_class(0,precision));
  reverse_copy(half_even.begin(), half_even.end(), even.begin());
  copy(half_even.begin(), half_even.end(), even.begin()+4+1);
    std::transform(odd.begin(), odd.begin() +4+1, odd.begin(), 
        [](mpf_class i) { 
          return mpf_class(-1.0*i,precision); 
        });
  for (mpf_class m : even) {
    cout << m  << "\n";
  }
  cout << endl;


}

void test_Gauss_Legendre_3()
{
  cout << "test gauss for legendre 3:\n";
  mpz_class highest_prec = 5;
  OrthogonalPolynomails* Pn = new LegendrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (auto gauss : gausstable) {
    cout << gauss;
  }
  cout << "3阶积分表完毕"<< endl;
}

void test_Gauss_Legendre_4()
{
  cout << "test gauss for legendre 4:\n";
  string dirpaths = "GaussPntDoc/Interval";
  mpz_class highest_prec = 7;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LegendrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (auto gauss : gausstable)
  {
    cout << gauss;
    writeGaussianInfo(dirpaths, gauss);
  }
  cout << "四阶积分表完毕\n";
}

void test_Legendre_100()
{
  cout << "test gauss for legendre 51:\n";
  string dirpaths = "GaussPntDoc/Interval";
  mpz_class highest_prec = 80;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LegendrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  //vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (int i = 1; i <= 41; ++i)
  {
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    writeGaussianInfo(dirpaths, gauss);
  }
  GaussianPoint1D gauss  = gtable.compute_gaussian_table(41);
  vector<mpf_class> coe(80,mpf_class(0,precision));
  coe.at(79) = 80;
  polynomial x_79(coe,79);
  vector<mpf_class> points = gauss.get_points();
  vector<mpf_class> weights = gauss.get_weights();
  mpf_class res(0,precision);
  for(int i = 0; i < points.size(); ++i)
  {
    res = res + weights.at(i)*x_79(points.at(i));
  }
  cout << "80*x^{79} 在[-1,1]上的积分为:\n";
  cout << fixed << setprecision(80) << res << endl;
  assert(abs(res) < 1e-80);
}
int main()
{
  //test_copy();
  //test_transform();
  //test_Gauss_Legendre_3();
  test_Gauss_Legendre_4();
  test_Legendre_100();
  return 0;
}
