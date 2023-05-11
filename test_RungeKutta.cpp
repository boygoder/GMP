#include "RungeKutta.h"
#include "PruferTransformer.h"
#include "OrthogonalPolynomials.h"
#include <iostream>
#include "math.h"
using namespace std;

void test_Prufer_for_Legendre2()
{
  cout << "test RungeKutta4 for Prufer Transform of Legendre2:\n";
  LegendrePolys Pn(mpz_class(2));
  polynomial p2 = Pn.getPm(2);
  vector<polynomial> coeff_ODE = Pn.coefficientsOfODE(2);
  PruferTransformer prufer2(p2,coeff_ODE);
  auto dx_dtheta = prufer2.get_dx_dtheta();
  vector<mpf_class> initial = {0,prufer2.theta_x_value(0)};
  RungeKutta4 RK4(initial,dx_dtheta);
  mpf_class h = -1e-3;
  mpf_class end_x = -pi/2.0;
  mpf_class end_y = RK4.compute(h, end_x);
  cout << "Legendre polynomial P2's smallest root>0 is :" << end_y <<"\n";
  cout << "should be about " << 1.0/sqrt(3) << endl;
  cout << "P2(" << end_y <<") = " << p2(end_y) << endl;
  cout << "P2(" << 1.0/sqrt(3) << ") = " << p2(1.0/sqrt(3)) << endl;
}

void test_Prufer_for_Legendre4()
{
  cout << "test RungeKutta4 for Prufer Transform of Legendre4:\n";
  LegendrePolys Pn(mpz_class(4));
  polynomial p4 = Pn.getPm(4);
  vector<polynomial> coeff_ODE = Pn.coefficientsOfODE(4);
  PruferTransformer prufer4(p4,coeff_ODE);
  auto dx_dtheta = prufer4.get_dx_dtheta();
  vector<mpf_class> initial = {0,prufer4.theta_x_value(0)};
  RungeKutta4 RK4(initial,dx_dtheta);
  mpf_class h = -1e-3;
  mpf_class end_x = -pi/2.0;
  mpf_class end_y = RK4.compute(h, end_x);
  cout << "Legendre polynomial P4's smallest root>0 is :" << end_y <<"\n";
  cout << "P4(" << end_y <<") = " << p4(end_y) << endl;
}

void test_set_initial()
{
  cout << "test RungeKutta4 for Legendre4 with set initial condtion:\n";
  LegendrePolys Pn(mpz_class(4));
  polynomial p4 = Pn.getPm(4);
  vector<polynomial> coeff_ODE = Pn.coefficientsOfODE(4);
  PruferTransformer prufer4(p4,coeff_ODE);
  auto dx_dtheta = prufer4.get_dx_dtheta();
  vector<mpf_class> initial = {0,prufer4.theta_x_value(0)};
  RungeKutta4 RK4(initial,dx_dtheta);
  mpf_class h = -1e-3;
  mpf_class end_x = -pi/2.0;
  mpf_class end_y = RK4.compute(h, end_x);
  cout << "Legendre polynomial P4's smallest root>0 is :" << end_y <<"\n";
  cout << "P4(" << end_y <<") = " << p4(end_y) << endl;

  vector<mpf_class> initial2 = {pi/2.0,end_y};
  RK4.set_initial_condition(initial2);
  mpf_class end_y2 = RK4.compute(h, end_x);
  cout << "Legendre polynomial P4's second root>0 is :" << end_y2 <<"\n";
  cout << "P4(" << end_y2 <<") = " << p4(end_y2) << endl;

}
int main()
{
  test_Prufer_for_Legendre2();
  test_Prufer_for_Legendre4();
  test_set_initial();
  return 0;
}
