#ifndef __POLYNOMIAL__  
#define __POLYNOMIAL__

//Thanks to Tingyuan Wang's polynomial class;
#include <iostream>
#include <vector>
#include <cmath>
#include "gmptools.h"
using namespace std;


class polynomial
{
private:
	vector<mpf_class> coeff = {mpf_class(0,precision)};
	mpz_class degree{0};
public:
	polynomial()
	{
		// vector<mpf_class> _coeff = {mpf_class(0,precision)};
		// coeff = _coeff;
	 //  degree = 0;
  }
	polynomial(vector<mpf_class> _coeff, mpz_class _degree)
    :coeff{_coeff},degree{_degree}
	{}
  polynomial(const polynomial& poly)
  {
    coeff = poly.coeff;
    degree = poly.degree;
  }
  ~polynomial() = default;
  void setCoeff(vector<mpf_class> _coeff)
  {
    coeff.clear();
    coeff.resize(_coeff.size());
    coeff = _coeff;
  }
  void setDegree(mpz_class _degree)
  {
    //if _degree is negative ,degree will be set to zero.\n";
    degree = max(_degree,mpz_class(0));
  }
  const vector<mpf_class>& getCoeff() const
  {
    return coeff;
  }
  const mpz_class& getDegree() const
  {
    return degree;
  }
  friend ostream& operator<<(ostream& out,const polynomial& p)
  {
    out << "degree of polynomial is:" << p.degree << endl;
    out << "coefficients of polynomials:\n";
    for(int i = 0; i < p.coeff.size(); ++i)
    {
      out << "degree " << i << " : " << p.coeff.at(i) << "\n";
    }
    return out;
  }
	friend polynomial operator+(polynomial a, polynomial b)
	{
		mpz_class a1(a.degree);
    vector<mpf_class> a2(a.coeff);
		mpz_class b1(b.degree);
		vector<mpf_class> b2(b.coeff);
		mpz_class n;
		if (a1 >= b1)
			n = a1;
		else
			n = b1;
		vector<mpf_class> x(n.get_ui()+1,mpf_class(0,precision));
		for (int i = 0; i <= n.get_ui(); i++)
		{
			if (i <= a1 && i <= b1)
				x[i] = a2[i] + b2[i];
			else if (i <= a1 && i > b1)
				x[i] = a2[i];
			else
				x[i] = b2[i];
		}
		polynomial c(x,n);
		return c;
	}
	friend polynomial operator*(mpf_class x, polynomial a)
	{
		mpz_class n(a.degree);
		vector<mpf_class> g(a.coeff);
		mpz_class m(n);
		vector<mpf_class> s(m.get_ui() + 1,mpf_class(0,precision));
		for (int i = 0; i <= m.get_ui(); i++)
			s[i] = x * g[i];
		polynomial c(s,m);
		return c;
	}
	friend polynomial operator-(polynomial a, polynomial b)
	{
		polynomial c = -1 * b;
		polynomial d = a + c;
		return d;
	}
	friend polynomial operator*(polynomial a, polynomial b)
	{
		mpz_class a1(a.degree);
		vector<mpf_class> a2(a.coeff);
		mpz_class b1(b.degree);
		vector<mpf_class> b2(b.coeff);
		mpz_class c = a1 + b1;
		vector<mpf_class> d(c.get_ui() + 1,mpf_class(0,precision));
		for (int i = 0; i <= c.get_ui(); i++)
		{
			d[i] = 0;
			for (int j = 0; j <= i; j++)
			{
				if (i - j <= b1 && j <= a1)
					d[i] += a2[j] * b2[i - j];
			}
		}
		polynomial e(d,c);
		return e;
	}
	// polynomial interpolate(mpf_class first, vector<mpf_class> root)
	// {
	// 	mpz_class n = root.size();
	// 	vector<mpf_class> a0 = {-root[0], 1.0};
	// 	vector<mpf_class> a1 = {-root[0], 1.0};
	// 	polynomial a;
	// 	polynomial b;
	// 	mpz_class p = 1;
	// 	a = polynomial(a0, p);
	// 	for (int i = 1; i < n.get_ui(); i++)
	// 	{
	// 		a1[0] = -root[i];
	// 		b = polynomial(a1, p);
	// 		a = a * b;
	// 	}
	// 	a = first * a;
	// 	return a;
	// }
	friend mpf_class pointValue(mpf_class x, polynomial a)
	{
		mpz_class a1 = a.degree;
		vector<mpf_class> a2(a.coeff);
		mpf_class c(0,precision);
		for (int i = 0; i <= a1.get_ui(); i++)
			c = c + a2[i] * mpf_class_pow_ui(x, i);
		return c;
	}
  mpf_class operator()(mpf_class x)
  {
		mpf_class c(0,precision);
		for (int i = 0; i <= degree.get_ui(); i++)
			c = c + coeff[i] * mpf_class_pow_ui(x, i);
		return c;
    //BUG: 秦九韶算法，存在精度问题
  //   vector<mpf_class> a = coeff;
  //   //reverse(a.begin(),a.end());
  //   vector<mpf_class> b(a.size(),mpf_class(0,precision));
  //   unsigned int order = degree.get_ui();
  //   b.at(order) = a.at(order);
		// for (int i = 1; i <= order; i++)
		// 	b.at(order-i) = b.at(order-(i-1))*x + a.at(order-i);
		// return b.at(0); 
  }
 mpf_class derivateValue(mpf_class x)
  {
    //TODO:秦九韶算法
    if (degree == 0) {
     return 0;
    }
    vector<mpf_class> a = coeff;
		vector<mpf_class> b(a.size(),mpf_class(0,precision));
    unsigned int order = degree.get_ui();
    b.at(order) = a.at(order);
		for (int i = 1; i <= order; i++)
			b.at(order-i) = b.at(order-(i-1))*x + a.at(order-i);
    vector<mpf_class> c(b.size()-1,mpf_class(0,precision));
    c.at(order-1) = b.at(order);
    for (int i = 1; i <= order-1; i++) {
      c.at(order-1-i) = c.at(order-i)*x + b.at(order-i);
    }
    return c.at(0);
  }
	friend polynomial derivate(polynomial a)
	{
		mpz_class a1(a.degree);
		vector<mpf_class> a2(a.coeff);
		mpz_class n;
		if (a1 == 0)
			n = 0;
		else
			n = a1 - 1;
		vector<mpf_class> x(n.get_ui() + 1,mpf_class(0,precision));
		if (a1 != 0) {
			for (int i = 0; i <= n.get_ui(); i++)
				x[i] = a2[i + 1] * (i + 1);
		}
		else {
			x[0] = 0;
    }
		polynomial e(x,n);
		return e;
	}
	// mpf_class integration(polynomial poly, mpf_class a, mpf_class b)
	// {
	// 	mpf_class key = 0;
	// 	vector<mpf_class> _coeff = poly.coeff;
	// 	mpz_class n = poly.degree;
	// 	for (mpz_class i = 0; i <= n; i++)
	// 	{
	// 		key += 1.0 * _coeff[i] * (pow(b, n + 1) - pow(a, n + 1)) / (i + 1);
	// 	}
	// 	return key;
	// }
};

#else
//do nothing
#endif
