#include <gmpxx.h>
#include <iostream>

using namespace std;
// g++ -std=c++17 -g test_gmp.cpp -o test_gmp -lgmp -lgmpxx
int main()
{
  mpz_class a, b, c;
  a = 1234;
  b = "-5678";
  c = a+b;
  cout << "sum is " << c << "\n";
  cout << "absolute value is " << abs(c) << "\n";
}
