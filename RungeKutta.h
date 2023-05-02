#ifndef __RUNGEKUTTA__
#define __RUNGEKUTTA__
#include <gmpxx.h>
#include <functional>
#include <vector>
#include <assert.h>
#include <iostream>
using namespace std;
class RungeKutta4
{
  private:
    using func = function<mpf_class(mpf_class,mpf_class)>;
    func derivate_f;
    vector<mpf_class> initial_condition;
    //mpf_class step_length;
    //vector<mpf_class> end_condition(2);
  public:
    RungeKutta4() = default;
    RungeKutta4(vector<mpf_class> initial,func f);
    ~RungeKutta4() = default;
    mpf_class compute(mpf_class h, mpf_class end_x);

};

RungeKutta4::RungeKutta4(vector<mpf_class> initial, func f)
{
  derivate_f = f;
  assert(initial.size() == 2);
  initial_condition = initial;
};

mpf_class RungeKutta4::compute(mpf_class h, mpf_class end_x)
{
  cout << "running RungeKutta4 method:\n";
  mpf_class x = initial_condition.at(0);
  mpf_class y = initial_condition.at(1);
  mpf_class length = end_x - x;
  assert(length*h > 0);
  mpf_class half_h = h/2.0;
  mpf_class one_sixth_h = h/6.0;
  mpf_class steps = length/h;
  cout << "initial condition is: (" << x <<"," << y << ")\n";
  cout << "step length is:" << h << ", num of steps is: " << steps <<"\n";
  cout << "end x is: "<< end_x <<"\n";
  mpf_class K1,K2,K3,K4;
  for(mpf_class i = 1; i <= steps; ++i)
  {
    K1 = derivate_f(x,y);
    K2 = derivate_f(x+half_h,y+half_h*K1);
    K3 = derivate_f(x+half_h,y+half_h*K2);
    K4 = derivate_f(x+h,y+h*K3);
    x = x + h;
    y = y + one_sixth_h*(K1+2*K2+2*K3+K4);
  }
  cout << "end y is: " << y << endl;
  return y;
};

#else
//do nothing
#endif
