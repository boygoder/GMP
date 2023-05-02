#include "NewtonMethod.h"
#include "OrthogonalPolynomials.h"
#include <iomanip>
#include "assert.h"
void test_Legendre2()
{
  cout << "test NewtonMethod for Legendre polynomials without talor:\n";
  LegendrePolys Pn(2);
  polynomial p2 = Pn.getPm(2);
  polynomial d_p2 = Pn.getDerivatePm(2);
  vector<polynomial> initial_function = {p2,d_p2}; 
  mpf_class start_x(0.57,precision);
  NewtonMethod n1(start_x,initial_function);
  mpf_class root16 = n1.compute_root(16);
  cout << "root of Legendre P2 with precision 16 is :\n" ;
  cout << "root16's precision is:" << root16.get_prec() << endl;
  cout << fixed << setprecision(80) << root16 << endl;
  cout << p2(root16) << endl;
  mpf_class root34 = n1.compute_root(34);
  cout << "root of Legendre P2 with precision 34 is :\n" ;
  cout << fixed << setprecision(80) <<root34 << endl;
  cout << p2(root34) << endl;
  //mpf_class root64 = n1.compute_root(64);
  //cout << "root of Legendre P2 with precision 64 is :\n" ;
  //cout << fixed << setprecision(80) <<root64 << endl;
  //cout << p2(root64) << endl;
  assert(abs(p2(root16)) < 1e-16);
  assert(abs(p2(root34)) < 1e-34);
  //assert(abs(p2(root64)) < 1e-64);
}
void test_Legendre2_taylor()
{
  cout << "test NewtonMethod for Legendre polynomials with taylor:" << endl;
  LegendrePolys Pn(2);
  polynomial p2 = Pn.getPm(2);
  polynomial d_p2 = Pn.getDerivatePm(2);
  mpf_class start_x(0.57,precision);
  vector<polynomial> initial_function = {p2,d_p2};
  vector<mpf_class> initial_value = {p2(start_x),d_p2(start_x)};
  vector<polynomial> coeff_OED = Pn.coefficientsOfODE(2);
  NewtonMethod n1(start_x,initial_function,coeff_OED);
  mpf_class root16 = n1.compute_root_with_taylor(30,16);
  cout << "root of Legendre P2 with precision 16 is :\n" ;
    cout << "root16's precision is:" << root16.get_prec() << endl;
  cout << fixed << setprecision(80) << root16 << endl;
  cout << p2(root16) << endl;
  mpf_class root34 = n1.compute_root_with_taylor(60,34);
  cout << "root of Legendre P2 with precision 34 is :\n" ;
  cout << fixed << setprecision(80) <<root34 << endl;
  cout << p2(root34) << endl;
  mpf_class root64 = n1.compute_root_with_taylor(120,64);
  cout << "precision of root64 is:" << root64.get_prec() << endl;
  cout << "root of Legendre P2 with precision 64 is :\n" ;
  cout << fixed << setprecision(80) <<root64 << endl;
  cout << p2(root64) << endl;
  mpf_class root80 = n1.compute_root_with_taylor(120,80);
  cout << "precision of root80 is:" << root80.get_prec() << endl;
  cout << "root of Legendre P2 with precision 80 is :\n" ;
  cout << fixed << setprecision(80) <<root80 << endl;
  cout << p2(root80) << endl;
  assert(abs(p2(root16)) < 1e-16);
  assert(abs(p2(root34)) < 1e-34);
  assert(abs(p2(root64)) < 1e-64);
  assert(abs(p2(root80)) < 1e-80);
}

int main()
{
  test_Legendre2();
  test_Legendre2_taylor();
  return 0;
};  
