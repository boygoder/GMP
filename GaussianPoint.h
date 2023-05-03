#ifndef __GAUSSIANPOINT__
#define __GAUSSIANPOINT__

#include "gmptools.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <filesystem>
#include <fstream>
#include <string>
using namespace std;
namespace fs = std::filesystem;
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

GaussianPoint1D::GaussianPoint1D(mpz_class _point_num, const vector<mpf_class>& _gauss_point, const vector<mpf_class>& _weight)
  :point_num{_point_num},gauss_point{_gauss_point},weight{_weight}
{
  assert(gauss_point.size() == weight.size());
};
vector<mpf_class>& GaussianPoint1D::get_points()
{
  return gauss_point;
};


vector<mpf_class>& GaussianPoint1D::get_weights()
{
  return weight;
};
void writeGaussianInfo(string dirpaths, const GaussianPoint1D& p)
{
    fs::path dirpath = dirpaths;
    if (!fs::exists(dirpath)) {
        fs::create_directory(dirpath);
        std::cout << "Directory created successfully.\n";
    }
    string filename = to_string(p.point_num.get_ui()) + ".txt"; 
    fs::path filepath = dirpath / filename;
    if (fs::exists(filepath)) {
        // 如果文件已经存在，则以覆盖模式打开
        std::ofstream ofs(filepath, ofstream::out | ofstream::trunc);
        if (ofs.is_open()) {
            // 文件成功打开，可以进行写操作
            ofs << p;
            ofs.close();
        } else {
            // 文件无法打开，输出错误信息并退出程序
            std::cerr << "Failed to open file " << filename <<".\n";
            exit(1);
        }
    } else {
        // 如果文件不存在，则创建文件
        std::ofstream ofs(filepath, ofstream::out);
        if (ofs.is_open()) {
            // 文件成功创建，可以进行写操作
            ofs << p;
            ofs.close();
        } else {
            // 文件无法创建，输出错误信息并退出程序
            std::cerr << "Failed to create file " << filename <<".\n";
            exit(1);
        }
    }


}
ostream& operator<<(ostream& out,const GaussianPoint1D& p)
{
  //out << "number of gauss points is: " << p.point_num << "\n";
  //out << "\t\tpoint" << "\t\t\t\t" << "weight" << "\n";
  out << p.point_num << "\n";
  vector<mpf_class> gauss_point = p.gauss_point;
  vector<mpf_class> weight = p.weight;
  for(int i = 0; i < gauss_point.size(); ++i)
  {
    out << fixed << setprecision(80);
    out << gauss_point.at(i) <<"\t";
    out << weight.at(i) << "\n";
  }
  out << endl;
  return out;
}
#else
// do nothing
#endif
