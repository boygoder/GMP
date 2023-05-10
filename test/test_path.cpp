#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
namespace fs = std::filesystem;   // 命名空间别名
using namespace std;
int main() {
    int a = 123;    // 指定数字 a 的值

    // 检查当前目录下是否有文件夹 'gauss'，如果没有则创建
    string dirpaths = "gauss";
    fs::path dirpath = dirpaths;
    if (!fs::exists(dirpath)) {
        fs::create_directory(dirpath);
        std::cout << "Directory created successfully.\n";
    }

    // 检查在 'gauss' 文件夹中是否有 'a.txt' 文件，如果没有则创建
    string filename = std::to_string(a) + ".txt"; 
    fs::path filepath = dirpath / filename;

    if (fs::exists(filepath)) {
        // 如果文件已经存在，则以覆盖模式打开
        std::ofstream ofs(filepath, std::ofstream::out | std::ofstream::trunc);
        if (ofs.is_open()) {
            // 文件成功打开，可以进行写操作
            ofs << "Hello, world!\n";
            ofs.close();
        } else {
            // 文件无法打开，输出错误信息并退出程序
            std::cerr << "Failed to open file " << filename <<".\n";
            exit(1);
        }
    } else {
        // 如果文件不存在，则创建文件
        std::ofstream ofs(filepath, std::ofstream::out);
        if (ofs.is_open()) {
            // 文件成功创建，可以进行写操作
            ofs << "Hello, world!\n";
            ofs.close();
        } else {
            // 文件无法创建，输出错误信息并退出程序
            std::cerr << "Failed to create file " << filename <<".\n";
            exit(1);
        }
    }




    std::cout << "Done.\n";
    return 0;
}

