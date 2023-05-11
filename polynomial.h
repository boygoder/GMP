#ifndef __POLYNOMIAL__  
#define __POLYNOMIAL__

//Thanks to Tingyuan Wang's polynomial class;
#include <vector>
#include "gmptools.h"
using namespace std;


class polynomial
{
private:
	vector<mpf_class> coeff = {mpf_class(0,precision)};
	mpz_class degree{0};
public:
	polynomial() = default;
	polynomial(vector<mpf_class> _coeff, mpz_class _degree);
  polynomial(const polynomial& poly);
  ~polynomial() = default;
  void setCoeff(vector<mpf_class> _coeff);
  void setDegree(mpz_class _degree);
  const vector<mpf_class>& getCoeff() const;
  const mpz_class& getDegree() const;
  mpf_class operator()(mpf_class x);
 	mpf_class derivateValue(mpf_class x);
	friend polynomial derivate(polynomial a);
	friend ostream& operator<<(ostream& out,const polynomial& p);
	friend polynomial operator+(polynomial a, polynomial b);
	friend polynomial operator*(mpf_class x, polynomial a);
	friend polynomial operator-(polynomial a, polynomial b);
	friend polynomial operator*(polynomial a, polynomial b);
	friend mpf_class pointValue(mpf_class x, polynomial a);
	// mpf_class integration(polynomial poly, mpf_class a, mpf_class b);
	// polynomial interpolate(mpf_class first, vector<mpf_class> root);
};

#else
//do nothing
#endif
