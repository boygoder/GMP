#include "PruferTransformer.h"
#include "OrthogonalPolynomials.h"
void test_Lengendre_2()
{
  LegendrePolys Pn(mpz_class(2));
  polynomial p2 = Pn.getPm(mpz_class(2));
  vector<polynomial> coeff_ODE = Pn.coefficientsOfODE(mpz_class(2));
  PruferTransformer prufer(p2,coeff_ODE);
  function<mpf_class(mpf_class)> theta_x = prufer.get_theta_x();
  cout << "for P2 polynomial of Legendre:\n";
  cout << "theta_x(0) = "<< theta_x(0) << "\n";
  function<mpf_class(mpf_class,mpf_class)> dx_dtheta = prufer.get_dx_dtheta();

  cout << "dx_dtheta(0,0) = " << dx_dtheta(0,0) << ", should be about " << -1.0/sqrt(6) << "\n";
  cout << "dx_dtheta(1,0) = " << dx_dtheta(1,0) << "\n";

  cout << "dx_dtheta(0,1) = " << dx_dtheta(0,1) << ",should be 0\n";
  cout << "dx_dtheta(1,1) = " << dx_dtheta(1,1) << "\n";
}

void test_Lengendre_3()
{
  LegendrePolys Pn(mpz_class(3));
  polynomial p3 = Pn.getPm(mpz_class(3));
  vector<polynomial> coeff_ODE = Pn.coefficientsOfODE(mpz_class(3));
  PruferTransformer prufer(p3,coeff_ODE);
  function<mpf_class(mpf_class)> theta_x = prufer.get_theta_x();
  cout << "for P3 polynomial of Legendre:\n";
  cout << "theta_x(0) = "<< theta_x(0) << "\n";
  function<mpf_class(mpf_class,mpf_class)> dx_dtheta = prufer.get_dx_dtheta();

  cout << "dx_dtheta(0,0) = " << dx_dtheta(0,0) << ",should be about " << -1/sqrt(12) <<"\n";
  cout << "dx_dtheta(1,0) = " << dx_dtheta(1,0) << "\n";

  cout << "dx_dtheta(0,1) = " << dx_dtheta(0,1) << "\n";
  cout << "dx_dtheta(1,1) = " << dx_dtheta(1,1) << "\n";
}


int main()
{
  test_Lengendre_2();
  test_Lengendre_3();
  return 0;
}
