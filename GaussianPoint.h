#ifndef __GAUSSIANPOINT__
#define __GAUSSIANPOINT__

#include "gmptools.h"
#include <vector>
#include <string>
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
    vector<mpf_class>& get_points();
    vector<mpf_class>& get_weights();
    friend void writeGaussianInfo(string dirpaths, const GaussianPoint1D& p);
    friend ostream& operator<<(ostream& out, const GaussianPoint1D& p); 
};

#else
// do nothing
#endif
