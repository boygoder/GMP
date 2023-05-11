#ifndef __PRUFERTRANSFORMER__
#define __PRUFERTRANSFORMER__

#include "polynomial.h"
#include <functional>
#include <vector>
using namespace std;

class PruferTransformer
{
  private:
    function<mpf_class(mpf_class)> theta_x;
    function<mpf_class(mpf_class,mpf_class)>  dx_dtheta;
  public:
    PruferTransformer() = default;
    PruferTransformer(polynomial& u, vector<polynomial>& coeff_ODE);
    ~PruferTransformer() = default;
    function<mpf_class(mpf_class)> get_theta_x();
    function<mpf_class(mpf_class,mpf_class)> get_dx_dtheta();
    mpf_class theta_x_value(mpf_class x);
    mpf_class dx_dtheta_value(mpf_class theta,mpf_class x);
};


#else
//do nothing
#endif
