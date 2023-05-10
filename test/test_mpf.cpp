#include <gmpxx.h>
#include <gmp.h>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
int main()
{
  mpf_class a(1.5,100);
  cout <<"value of a:" << a << "\n";
  cout << "precision of a is:" << a.get_prec() << "\n";
  mpf_class b[1];
  for(int i = 0; i < 1; ++i) {
    b[i].set_prec(a.get_prec());
    cout << "precision of b is:" << b[i].get_prec() << "\n";
    b[i] = a;
    cout << "precision of b is:" << b[i].get_prec() << "\n";
  }
  vector<mpf_class> vf(1,mpf_class(3,64));
  cout << "precision of vf[0] is:" << vf[0].get_prec() << " ,value is :" << vf[0] <<"\n";
  using Real = mpf_class;
  vector<Real> vf2(10,mpf_class(2));
  cout << vf2.size() << "\n";
  cout << vf2[0] << "\n";
  unsigned int ex = 3;
  cout << vf2[0].get_mpf_t() << "\n";
  mpz_class base=2;
  mpz_t base_t;
  mpz_init_set(base_t,base.get_mpz_t());
  mpz_t result_t;
  mpz_init(result_t);
  mpz_pow_ui(result_t,base_t,10);
  cout << base_t <<" " << result_t <<"\n";
  cout << endl;
  mpz_clear(base_t);
  mpz_clear(result_t);
}
