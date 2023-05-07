#include "gmptools.h"
#include <iostream>
#include "OrthogonalPolynomials.h"
#include <iomanip>

void test_prec()
{
  mpf_class a1(0,256);
  mpf_class a2(0);
  mpf_class a3(0,128);
  cout << "a1's precision:" << a1.get_prec() << endl;
  cout << "a2's precision:" << a2.get_prec() << endl;
  cout << "a3's precision:" << a3.get_prec() << endl;
  cout << "test assignment:\n";
  mpf_class b = a1;
  cout << "assigment construct:" << b.get_prec() << endl;
  a1 = a2;
  cout << "after a1=a2,a1's precision is:" << a1.get_prec()<<endl;
  a2 = a3;
  cout << "after a2=a3,a2's precision is:" << a2.get_prec()<<endl;
  cout << "test construct:\n";
  mpf_class a4(0,128);
  mpf_class a5(0,32);

  cout << "a4's precision:" << a4.get_prec() << endl;
  cout << "a5's precision:" << a5.get_prec() << endl;
  mpf_class a6(a5);
  cout << "after a6(a5),a6's precision is:" << a6.get_prec() << endl;
  mpf_class a7(a4);
  cout << "after a7(a4),a7's precision is:" << a7.get_prec() << endl;

}

void test_vector()
{
  cout << "test vector precision:\n";
  vector<mpf_class> a1(1);
  cout << "vector default precision:" << a1[0].get_prec() << endl;
  vector<mpf_class> a2 = {mpf_class(0,256)};
  a2.push_back(1);
  cout << "first is 256,the second is:" << a2[1].get_prec() << endl;
  vector<mpf_class> a = {mpf_class(0,128),mpf_class(1,128)};
  cout << "after b(a), b[0]'s precision is:";
  vector<mpf_class> b(a);
  cout << b[0].get_prec() << endl;

  vector<mpf_class> c = a;
  cout << "after c=a, c[0]'s precision is:" << c[0].get_prec() << endl;

  vector<mpf_class> d(a.begin(),a.end());
  cout << "after d(a.begin(),a.end()),d[0]'s precision is:" << d[0].get_prec() << endl;
}

void test_add()
{
  mpf_class a(1,256);
  mpf_class b(0.5,128);
  a = a+b;
  cout << "after a=a+b, a's precision is:" << a.get_prec()<< endl;
}

void test_polynomial()
{
  LegendrePolys Pn(2);
  polynomial p2 = Pn.getPm(2);
  mpf_class x(0.5,precision);
  mpf_class res = p2(x);
  cout << "res's precision is:" << res.get_prec() << endl;
  mpf_class x2(0.5);
  mpf_class res2 = p2(x2);
  cout << " when x's precision is 64,res's precision is:"<< res2.get_prec() << endl;
  polynomial p1 = Pn.getPm(1);
  polynomial p1_2 = p1 + p2;
  mpf_class res3 = p1_2(x2); 
  cout << "p1+p2 res's precision is:" << res3.get_prec() << endl;
}

void test_mpz_div()
{
  mpz_class prec = 100;
  mpz_class order = (prec+1)/2;
  cout << "(100+1)/2 = "<< order << endl;
  mpf_class order2 = (prec+1)/2.0;
  cout << "mpf (100+1)/2 = " << order2 << endl;
  mpf_class order3 = ceil((prec.get_d()+1)/2.0);
  cout << "ceil get_d (100+1)/2 = "<< order3 << endl;
}

void test_vector_construct() {
  cout << "test vector construct:\n";
  int len = 3;
  cout << "vector v{int,mpf}:\n";
  vector<mpf_class> v1{len,mpf_class(0,precision)};
  for(auto i : v1){
    cout << i << endl;
  }
  cout << "vector v(int,mpf):\n";
  vector<mpf_class> v2(len,mpf_class(0,precision));
  for (auto i : v2){
    cout << i << endl;
  }

}

void test_pi()
{
  cout << "now pi is :\n";
  cout << fixed << setprecision(128) << pi << "\n";
  cout << "new pi is :\n";
  mpf_class new_pi = compute_pi(128);
  cout << new_pi << endl;
}
int main()
{
  test_vector();
  test_prec();
  test_add();
  test_polynomial();
  test_mpz_div();
  test_vector_construct();
  test_pi();
  return 0;
}
