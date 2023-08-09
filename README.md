# GMP

## 概述

利用gmp库提供的高精度浮点数，通过求解勒让德多项式、拉盖尔多项式、赫米特多项式的根，得到高精度的高斯积分表，对于数值计算具有重要的意义。


## gmp安装

[官网](https://gmplib.org/)，默认安装在/usr/local/

```bash
# 安装前准备：（安装m4）
sudo apt install m4
wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz
tar xvJf ***.tar.xz
cd gmp-6.2.1
sudo ./configure --enable-cxx=detect
make && make check &&sudo make install
```

引用头文件

```C++
//C头文件
#include "gmp.h"
//编译
gcc mycprog.c  -o myprog  -lgmp
//C++头文件
#include <gmpxx.h>
//编译
g++ mycxxprog.cc  -o mycxxprog -lgmpxx -lgmp
```


## 数值结果

### 勒让德多项式

对于偶数阶多项式，$x^{n}$在$[-1,1]$上的积分为$\frac{2}{n+1}$。

表格中是数值积分值与理论积分值的误差的绝对值。

| 测试多项式\浮点数精度 |  $10^{-16}$  |  $10^{-34}$   |  $10^{-80}$    |
| :--: | :--: | :--: | :--: |
|    $x^{20}$    | 9.9e-20 | 4.6e-41 | 1.7e-117 |
|    $x^{50}$    | 3.8e-19 | 1.3e-37 | 5.7e-108 |
|    $x^{100}$    | 1.4e-19 | 1.7e-38 | 4.4e-105 |

### 赫米特多项式

此处理论值采用从wolframalpha得到的高精度值。

| 测试多项式\浮点数精度 |  $10^{-16}$  |  $10^{-34}$   |  $10^{-80}$    |
| :--: | :--: | :--: | :--: |
|    $x^{20}$    | 3.3e-21 | 6.4e-48 | 2.9e-101 |
|    $x^{50}$    | 1.5e-29 | 4.7e-34 | 2.9e-83 |
|    $x^{100}$    | 2.2e-16 | 6.4e-44 | 1.6e-82 |

