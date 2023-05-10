#include <gmpxx.h>
#include <iostream>
using namespace std;
int main()
{
  mpf_class a(1.5,100);
  mpf_t neg_a_t;
  mpf_init(neg_a_t);
  mpf_neg(a.get_mpf_t(),a.get_mpf_t());
  // mpf_neg(neg_a_t,a.get_mpf_t());
  mpf_class neg_a(neg_a_t,a.get_prec());
  cout <<"neg_a:" << neg_a <<", " << neg_a.get_prec() << endl;
  cout <<"a:" << a <<"," <<  a.get_prec() << endl;
}
