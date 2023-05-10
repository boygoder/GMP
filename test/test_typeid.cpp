#include <iostream>
#include <typeinfo>
#include <assert.h>
#include <memory>
using namespace std;


class A
{
  public:
    virtual void name() = 0;
};

class B: public A
{
  public:
    virtual void name() {
      cout << "B" << endl;
    }
};

class C :public A
{
  public:
    virtual void name() {
      cout << "C" << endl;
    }
};
int main()
{
	A* vc = new C();
  //必须存在虚函数，才能判断为派生类。
	assert(typeid(*vc) == typeid(C));
  assert(typeid(*vc) != typeid(B));
	cout << typeid(*vc).name() << endl;
	A* vb = new B();
	assert(typeid(*vb) == typeid(B));
	cout << typeid(*vb).name() << endl;
	unique_ptr<A> uvb(vb);
  cout << typeid(*uvb.get()).name() << endl;
  return 0;
}
