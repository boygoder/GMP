#include "GaussIntegralTableGenerator.h"
#include <iostream>
#include "math.h"
#include "assert.h"
#include <iomanip>
using namespace std;

void Gauss_Legendre_Table(mpz_class highest_prec)
{
  cout << "generate gauss integral table via Legendre polynomials……:\n";
  string dirpaths = "GaussPntDoc/Interval/Legendre";
  // 对应的多项式阶数，n个点可达到2n-1精度。
  mpz_class highest_order = ceil((highest_prec.get_d()+1)/2.0);
  //选择正交多项式
  OrthogonalPolynomials* Pn = new LegendrePolys();
  //建立高斯积分表生成类
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  //设置泰勒展开项数和牛顿迭代法精度
  gtable.set_taylor_items(120);
  gtable.set_newton_eps(80);
  //从1阶到highest_order阶的积分信息生成，并保存到文件。
  for (int i = 1; i <= highest_order.get_ui(); ++i)
  {
    cout << "compute order " << i << endl;
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    writeGaussianInfo(dirpaths, gauss);
  }

}


void Gauss_Hermite_Table(mpz_class highest_prec)
{
  cout << "generate gauss integral table via Hermite polynomials……:\n";
  string dirpaths = "GaussPntDoc/Interval/Hermite";
  mpz_class highest_order = ceil((highest_prec.get_d()+1)/2.0);
  OrthogonalPolynomials* Pn = new HermitePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  gtable.set_taylor_items(120);
  gtable.set_newton_eps(80);
  for (int i = 1; i <= highest_order.get_ui(); ++i)
  {
    cout << "compute order " << i << endl;
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    writeGaussianInfo(dirpaths, gauss);
  }
}

void Gauss_Lagueree_Table(mpz_class highest_prec)
{
  cout << "generate gauss integral table via Lagueree polynomials……:\n";
  string dirpaths = "GaussPntDoc/Interval/Lagueree";
  mpz_class highest_order = ceil((highest_prec.get_d()+1)/2.0);
  OrthogonalPolynomials* Pn = new LaguerrePolys();
  GaussIntegralTableGenerator gtable(highest_prec,Pn);
  gtable.set_taylor_items(120);
  gtable.set_newton_eps(80);
  for (int i = 1; i <= highest_order.get_ui(); ++i)
  {
    cout << "compute order " << i << endl;
    GaussianPoint1D gauss = gtable.compute_gaussian_table(i);
    cout << gauss;
    writeGaussianInfo(dirpaths, gauss);
   }
}
int main()
{
  int prec = 1;
  cout << "please input the prec you want:\n";
  cin >> prec;
  Gauss_Legendre_Table(prec);
  Gauss_Hermite_Table(prec);
  Gauss_Lagueree_Table(prec);
  return 0;

}
