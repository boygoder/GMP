#ifndef __GAUSSIANPOINT__
#define __GAUSSIANPOINT__

#include "gmptools.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <assert.h>
using namespace std;

class GaussianPoint1D
{
  private:
    vector<mpf_class> gauss_point;
    vector<mpf_class> weight;
    mpz_class point_num;
  public:
    GaussianPoint1D() = default;
    GaussianPoint1D(mpz_class _point_num,const vector<mpf_class>& _gauss_point, const vector<mpf_class>& _weight);
    ~GaussianPoint1D() = default;
    void writeGaussianInfo(string filename);
    friend ostream& operator<<(ostream& out, const GaussianPoint1D& p); 
};

GaussianPoint1D::GaussianPoint1D(mpz_class _point_num, const vector<mpf_class>& _gauss_point, const vector<mpf_class>& _weight)
  :point_num{_point_num},gauss_point{_gauss_point},weight{_weight}
{
  assert(gauss_point.size() == weight.size());
};
ostream& operator<<(ostream& out,const GaussianPoint1D& p)
{
  out << "number of gauss points is: " << p.point_num << "\n";
  out << "\t\t\tpoint" << "\t\t\t\t\t\t" << "weight" << "\n";
  vector<mpf_class> gauss_point = p.gauss_point;
  vector<mpf_class> weight = p.weight;
  for(int i = 0; i < gauss_point.size(); ++i)
  {
    out << fixed << setprecision(32);
    out << gauss_point.at(i) <<"\t";
    out << weight.at(i) << "\n";
  }
  out << endl;
  return out;
}
#else
// do nothing
#endif
