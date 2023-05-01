#include "OrthogonalPolynomials.h"
#include "polynomial.h"
#include <gmpxx.h>


void test_Legendre()
{
  mpz_class hdegree(4);
  vector<polynomial> theory_p(5);
  for (int i = 0; i <= hdegree.get_ui(); ++i)
  {
    theory_p.at(i).setDegree(mpz_class(i));
  }
  vector<mpf_class> coeff_p0 = {1};
  theory_p.at(0).setCoeff(coeff_p0);
  vector<mpf_class> coeff_p1 = {0,1};
  theory_p.at(1).setCoeff(coeff_p1);
  vector<mpf_class> coeff_p2 = {-0.5,0,1.5};
  theory_p.at(2).setCoeff(coeff_p2);
  vector<mpf_class> coeff_p3 = {0,-1.5,0,2.5};
  theory_p.at(3).setCoeff(coeff_p3);
  vector<mpf_class> coeff_p4 = {3.0/8,0,-30.0/8,0,35.0/8};
  theory_p.at(4).setCoeff(coeff_p4);

  vector<polynomial> theory_p_d(5);
  for (int i = 0; i <= hdegree.get_ui(); ++i)
  {
    theory_p_d.at(i).setDegree(mpz_class(i-1));
  }
  vector<mpf_class> coeff_p0_d = {0};
  theory_p_d.at(0).setCoeff(coeff_p0_d);
  vector<mpf_class> coeff_p1_d = {1};
  theory_p_d.at(1).setCoeff(coeff_p1_d);
  vector<mpf_class> coeff_p2_d = {0,3};
  theory_p_d.at(2).setCoeff(coeff_p2_d);
  vector<mpf_class> coeff_p3_d = {-1.5,0,7.5};
  theory_p_d.at(3).setCoeff(coeff_p3_d);
  vector<mpf_class> coeff_p4_d = {0,-15.0/2,0,35.0/2};
  theory_p_d.at(4).setCoeff(coeff_p4_d);



  LegendrePolys Pn(hdegree);
  polynomial p0 = Pn.getPm(0);
  cout << "compute p0:" << p0 << endl;
  polynomial p1 = Pn.getPm(1);
  cout << "compute p1:" << p1<< endl;
  polynomial p2= Pn.getPm(2);
  cout << "compute p2:" << p2 << endl;
  cout << "theory p2:" << theory_p.at(2) << endl;
  polynomial p3 = Pn.getPm(3);
  cout << "compute p3:" << p3 << endl;
  cout << "theory p3:" << theory_p.at(3) << endl;
  polynomial p4 = Pn.getPm(4);
  cout << "compute p4:" << p4 << endl;
  cout << "theory p4:" << theory_p.at(4) << endl;

  polynomial p0_d = Pn.getDerivatePm(0);
  cout << "compute p0_d:" << p0_d << endl;
  polynomial p1_d = Pn.getDerivatePm(1);
  cout << "compute p1_d:" << p1_d << endl;
  polynomial p2_d = Pn.getDerivatePm(2);
  cout << "compute p2_d:" << p2_d << endl;
  cout << "theory p2_d:" << theory_p_d.at(2) << endl;
  polynomial p3_d = Pn.getDerivatePm(3);
  cout << "compute p3_d:" << p3_d << endl;
  cout << "theory p3_d:" << theory_p_d.at(3) << endl;
  polynomial p4_d = Pn.getDerivatePm(4);
  cout << "compute p4_d:" << p4_d << endl;
  cout << "theory p4_d:" << theory_p_d.at(4) << endl;
}

int main()
{
  test_Legendre();
  return 0;
}
