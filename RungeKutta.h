#ifndef __RUNGEKUTTA__
#define __RUNGEKUTTA__
#include "gmptools.h"
#include <functional>
#include <vector>
using namespace std;
class RungeKutta4
{
  private:
    using func = function<mpf_class(mpf_class,mpf_class)>;
    func derivate_f;
    vector<mpf_class> initial_condition;
  public:
    RungeKutta4();
    RungeKutta4(vector<mpf_class> initial,func f);
    ~RungeKutta4() = default;
    mpf_class compute(mpf_class h, mpf_class end_x);
    void set_initial_condition(vector<mpf_class> initial);
};


#else
//do nothing
#endif
