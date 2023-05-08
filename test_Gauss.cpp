#include <gmp.h>
#include <iomanip>
#include <vector>
#include "GaussianPoint.h"
#include "OrthogonalPolynomials.h"
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
  string dirpaths = "GaussPntDoc/Interval/Legendre";
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

void test_Legendre_integral()
{
  cout << "test gauss integral for legendre :\n";
  mpz_class highest_prec = 100;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LegendrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  GaussianPoint1D gauss  = gtable.compute_gaussian_table(51);
  vector<mpf_class> coe(101,mpf_class(0,precision));
  coe.at(100) = 1;
  polynomial x_100(coe,100);
  vector<mpf_class> points = gauss.get_points();
  vector<mpf_class> weights = gauss.get_weights();
  mpf_class res(0,precision);
  for(int i = 0; i < points.size(); ++i)
  {
    res = res + weights.at(i)*x_100(points.at(i));
  }
  cout << "x^{100} 在[-1,1]上的积分为:\n";
  cout << fixed << setprecision(256) << res << endl;
  mpf_class theory = mpf_class(2.0,precision)/mpf_class(101,precision);
  cout << "理论积分为:\n";
  cout << fixed << setprecision(256) << theory << endl;
  cout << "误差为:\n";
  cout << fixed << setprecision(256) << res-theory << endl;  
  assert(abs(res-theory) < 1e-80);
}

void test_Legendre_100()
{
  cout << "test gauss for legendre 51:\n";
  string dirpaths = "GaussPntDoc/Interval/Legendre";
  mpz_class highest_prec = 100;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LegendrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  // vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (int i = 1; i <= 51; ++i)
  {
    cout << "compute order " << i << endl;
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    writeGaussianInfo(dirpaths, gauss);
  }

}

void test_Legendre_200()
{
  cout << "test gauss for Legendre 100:\n";
  string dirpaths = "GaussPntDoc/Interval/Legendre";
  mpz_class highest_prec = 199;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LegendrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  // vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (int i = 100; i <= 100; ++i)
  {
    cout << "compute order " << i << endl;
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    // writeGaussianInfo(dirpaths, gauss);
  }
}

void test_Hermite_200()
{
  cout << "test gauss for Hermite 100:\n";
  string dirpaths = "GaussPntDoc/Interval/Hermite";
  mpz_class highest_prec = 199;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new HermitePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  gtable.set_taylor_items(120);
  gtable.set_newton_eps(80);
  // vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (int i = 100; i <= 100; ++i)
  {
    cout << "compute order " << i << endl;
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    // writeGaussianInfo(dirpaths, gauss);
  }
}


void test_Hermite_integral()
{
 cout << "test gauss integral for Hermite :\n";
  mpz_class highest_prec = 100;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new HermitePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  GaussianPoint1D gauss  = gtable.compute_gaussian_table(51);
  vector<mpf_class> coe(101,mpf_class(0,precision));
  coe.at(100) = 1;
  polynomial x_100(coe,100);
  vector<mpf_class> points = gauss.get_points();
  vector<mpf_class> weights = gauss.get_weights();
  mpf_class res(0,precision);
  for(int i = 0; i < points.size(); ++i)
  {
    res = res + weights.at(i)*x_100(points.at(i));
  }
  cout << "x^{100} e^{-x^2} 在[-infinity,infinity]上的积分为:\n";
  cout << fixed << setprecision(256) << res << endl;
  mpf_t theory_t;
  mpf_init2(theory_t,precision);
  mpf_set_str(theory_t,"4.2904629123519598109157551960589376738242902258244938240430485310826851318323799037707311554284993957036806886756317727235137848794812440419209013632578077067211610203085038864055272776263386097646164e63",10);
  mpf_class theory(theory_t);
  cout << "理论积分为:\n";
  cout << fixed << setprecision(256) << theory << endl;
 cout << "误差为:\n";
  cout << fixed << setprecision(256) << res-theory << endl;  
  assert(abs(res-theory) < 1e-80);
  mpf_clear(theory_t);
}

void test_Lagueree_200()
{
  cout << "test gauss for Lagueree 100:\n";
  string dirpaths = "GaussPntDoc/Interval/Lagueree";
  mpz_class highest_prec = 199;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LaguerrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  gtable.set_taylor_items(120);
  gtable.set_newton_eps(80);
  // vector<GaussianPoint1D> gausstable = gtable.compute_gaussian_table();
  for (int i = 100; i <= 100; ++i)
  {
    cout << "compute order " << i << endl;
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    // writeGaussianInfo(dirpaths, gauss);
   }
}
void test_Lagueree_integral()
{
 cout << "test gauss integral for Lagueree :\n";
  mpz_class highest_prec = 100;
  // should call set_highest_degree
  OrthogonalPolynomails* Pn = new LaguerrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  GaussianPoint1D gauss  = gtable.compute_gaussian_table(31);
  vector<mpf_class> coe(61,mpf_class(0,precision));
  coe.at(60) = 1;
  polynomial x_100(coe,60);
  vector<mpf_class> points = gauss.get_points();
  vector<mpf_class> weights = gauss.get_weights();
  mpf_class res(0,2*precision);
  for(int i = 0; i < points.size(); ++i)
  {
    mpf_class ux(0,2*precision);
    ux = x_100(points.at(i));
    mpf_class wx(0,2*precision);
    wx = weights.at(i);
    res = res + wx*ux;
  }
  cout << "x^{60} e^{-x} 在[0,infinity]上的积分为:\n";
  cout << fixed << setprecision(256) << res << endl;
  mpf_t theory_t;
  mpf_init2(theory_t,precision);
  mpf_set_str(theory_t,"8320987112741390144276341183223364380754172606361245952449277696409600000000000000",10);
/*   mpf_set_str(theory_t,"9.3326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864e157",10); */
  mpf_class theory(theory_t,precision);
  cout << "理论积分为:\n";
  cout << fixed << setprecision(256) << theory << endl;
 cout << "误差为:\n";
  cout << fixed << setprecision(256) << res-theory << endl;  
  assert(abs((res-theory)/theory) < 1e-80);
  mpf_clear(theory_t);
}
int main()
{
  //test_copy();
  //test_transform();
  //test_Gauss_Legendre_3();
  //test_Gauss_Legendre_4();
  // test_Legendre_integral();
  //test_Legendre_100();
  // test_Legendre_200();
  // test_Hermite_200();
  // test_Hermite_integral();
  // test_Lagueree_200();
  test_Lagueree_integral();
  return 0;

}
