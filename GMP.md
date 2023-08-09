# GMP

## 安装

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

## 算法

### 主题

本文的主题是观察到经典特殊函数的根，可以通过利用这些函数满足如下形式的常微分方程来求解。
$$
p(x) \dfrac{d^2 u}{dx^2}(x) + q(x) \dfrac{du}{dx}(x) + r(x)u(x) = 0
\tag{1}
\label{eq1}
$$
其中$p(x),q(x),r(x)$是二阶多项式。

通过对$\eqref{eq1}$式进行微分，可以得到$u(x)$的高阶导数值，
$$
\begin{aligned}
p u^{(k+2)} 
&= -(kp' + q) u^{(k+1)}  \\
&- (\dfrac{k(k-1)}{2}p'' + kq' + r) u^{(k)} \\
&- (\dfrac{k(k-1)}{2}q'' + kr') u^{(k-1)} \\
&-\dfrac{k(k-1)}{2} r'' u^{(k-2)}
\end{aligned}
\tag{10}
\label{eq10}
$$
**注意**：此处我们所得到的递推式式针对于函数整体的连续导数，而在后面的牛顿法中，我们所需要的初始条件含有泰勒展开的单一项，需要做简单的处理。

本文的数值程序分为两步去求解多项式的每一个根，但是方程的第一个根需要其他手段来确定。在第一步中，通过辅助微分方程计算出一个根的逼近，使用的方法是Runge-Kutta法。辅助方程通过Prufer变换得到，初值由上一个计算得到的根确定。在第二步中，通过牛顿法和原多项式的泰勒级数展开格式，提高根的精度。

### Prufer变换

给定$\eqref{eq1}$式和可微正函数$\gamma : \mathbb{R} \rightarrow \mathbb{R}$，我们定义函数$\theta : \mathbb{R} \rightarrow \mathbb{R}$作为下列微分方程的解：
$$
\theta' = -\dfrac{r}{p} \sin^2(\theta) - \dfrac{r}{\gamma} \cos^2(\theta) - (\dfrac{\gamma'}{\gamma} + \dfrac{q-p'}{p})\dfrac{sin(2\theta)}{2}
\tag{2}
\label{eq2}
$$
其中$\theta,\gamma,p,q,r$是关于$x$的函数。

我们得到$\theta(x)$的表达式为:
$$
\theta(x) = \arctan(\dfrac{1}{\gamma(x)}  \dfrac{p(x) u'(x)}{u(x)}) + k\pi
\tag{3}
\label{eq3}
$$
其中$u(x)$即为满足$\eqref{eq1}$式的多项式，k是任意整数。

观察到对于任意满足$u(\tilde{x}) = 0$的$\tilde{x}$，$\theta(\tilde{x}) = (n+1/2)\pi$,n为整数。同样地，对于任意满足$u'(\tilde{x}) = 0$的$\tilde{x}$，$\theta(\tilde{x}) = n\pi$,n为整数。

为了让$\theta(x)$关于x的变化更加光滑和便于计算，选择$\gamma(x) = \sqrt{rp} $,于是$\eqref{eq2}$式和$\eqref{eq3}$式写为：
$$
\theta(x) = \arctan(\dfrac{1}{\sqrt{r(x)p(x)}}  \dfrac{p(x) u'(x)}{u(x)}) 
\tag{4}
\label{eq4}
$$

$$
\theta' = -\sqrt{\dfrac{r}{p}} - \dfrac{r'p - p'r + 2rq}{2rp} \dfrac{\sin(2\theta)}{2} 
\tag{5}
\label{eq5}
$$

对于特殊函数，都满足$\dfrac{dx}{d\theta} < 0$,$\eqref{eq5}$可以重写为如下方程，将x看作是关于$\theta$的函数，

$$
\dfrac{dx}{d\theta} = -(\sqrt{\dfrac{r}{p}} + \dfrac{r'p - p'r + 2rq}{2rp} \dfrac{\sin(2\theta)}{2} )^{-1}
\tag{6}
\label{eq6}
$$




### 正交多项式

**定义**：如果对于给定的区域(a,b)和正权重函数$\omega$,若对任意的k阶多项式$p_k$总能满足：
$$
\int_{a}^{b} p_{i}(x)p_{j}(x) \omega(x) dx 
= \left\{
\begin{aligned}
&0，i \neq  j \\  
&a_i, i = j 
\end{aligned}
\right.
\tag{7}
\label{eq7}
$$
则称$p_i$是区间（a,b）上关于权函数$\omega$构成正交族的一组多项式。

勒让德、赫米特和拉盖尔多项式都有着如下的循环关系：
$$
a_m p_{m+1}(x) = (b_m + c_m x) p_{m}(x) - d_m p_{m-1}(x)
\tag{8}
\label{eq8}
$$
对其求导可以得到导数的循环关系：
$$
a_m p'_{m+1}(x) = (b_{m} + c_{m}x)p'_{m}(x) - d_m p'_{m-1}(x) + c_{m} p_{m}(x)
\tag{9}
\label{eq9}
$$
其中$a_m,b_m,c_m,d_m$见下表，$P_n$表示勒让德多项式，$H_n$表示赫米特多项式，$L_n$表示拉盖尔多项式。

![image-20230506164917839](./assets/image-20230506164917839.png)

定义域$[a,b]$和权重$\omega$、对应的常微分方程系数$p,q,r$如下：

![image-20230506165210952](./assets/image-20230506165210952.png)

- 勒让德多项式：区间为[-1,1]，权函数$w=1$
  $$
  P_{0}(x ) = 1,P_{n}(x)  = \dfrac{1}{2^n n!} \dfrac{d^{n}}{dx^{n}} \{(x^2 - 1)^{n}\}(n=1,2,\cdots)
  \nonumber
  $$
  递推公式为：$(n+1)P_{n+1}(x) = (2n+1)xP_{n}(x)-nP_{n-1}(x)$,$P_0(x) =1,P_{1}(x) = x$.
  
  n为偶数时，Pn是偶函数，x=0是一个极值点，利用算法2求第一个根。
  
  n为奇数时，Pn是奇函数，x=0是第一个根。
  
- 赫米特多项式：区间为$[-\infty,\infty]$,权函数为$e^{-x^2}$,
  $$
  H_{0}(x)= 1,H_{n}(x) = (-1)^{n} e^{x^2} \dfrac{d^n}{dx^n}(e^{-x^2})
  \nonumber
  $$
  
  
  递推公式为$H_{0}(x) =1,H_{1}(x) = 2x,H_{n+1} (x) = 2xH_{n}(x) -2n H_{n-1}(x)$.
  
  n为偶数时，Pn是偶函数，x=0是一个极值点，利用算法2求第一个根。
  
  n为奇数时，Pn是奇函数，x=0是第一个根。
  
- 拉盖尔多项式：区间为$[0,\infty]$,权函数为$e^{-x}$,

  

$$
L_{0}(x)=1,L_{n}(x) = \dfrac{e^x}{n!} \dfrac{d^n}{dx^n}(x^n e^{-x})
\nonumber
$$
递推式为$L_{0}(x)=1,L_{1}(x) = 1-x,(n+1)L_{n+1}(x) = (2n+1-x)L_{n}(x) - nL_{n-1}(x)$.

根据递归式$x L_{n}'(x) = nL_{n}(x) - n L_{n-1}(x) = (x-n-1)L_{n}(x) +(n+1)L_{n+1}(x)$,当x是$L_{n}(x)$的一个根时，$xL_{n}'(x) = -nL_{n-1}(x) = (n+1)L_{n+1}(x)$。

拉盖尔多项式的权重$w_{n}^{k} = \dfrac{1}{x_{k} L_{n}'^{2}(x_{k})} = \dfrac{x_{k}}{(n+1)^2 L_{n+1}^{2}(x_k)} = \dfrac{x_{k}}{n^2 L_{n-1}^{2}(x_{k})}$。

拉盖尔多项式不具有奇偶性，但是n阶拉盖尔多项式的最小根满足$x_{1}^{n} > \dfrac{2}{4n+2}$,利用算法3求第一个根。

**如何赫米特多项式和拉盖尔多项式可能会溢出，可以使用正则化方式。**



给定区间(a,b)上，权重为$\omega$的一个正交族的n阶多项式的n个根$x_1,\cdots,x_n$,存在实数$w_{1}^{n},\cdots,w_{n}^{n}$，使得$\int_{a}^{b} f(x)\omega(x)dx \approx \sum_{k=1}^{n} \omega_{n}^{k} f(x_{k})$.对于阶数低于2n-1的多项式f是精确的。

### 龙格库塔法

![image-20230430165022051](./assets/image-20230430165022051.png)

### 牛顿法

设已知方程$f(x) = 0$有近似根$x_k$(假定$f'(x_k) \neq 0$),将函数$f(x)$在点$x_k$展开，有
$$
f(x) \approx f(x_k) + f'(x_k)(x-x_k)
$$
于是方程$f(x) = 0$可以近似地表示为
$$
f(x_k) + f'(x_k)(x-x_k) = 0
$$
这是个线性方程，记其根为$x_{k+1}$，则$x_{k+1}$的计算公式为
$$
x_{k+1} = x_{k} - \dfrac{f(x_k)}{f'(x_k)}
$$
这就是牛顿法，具有二阶收敛性。

对于正交多项式方程$u(x)=0$，在计算出$x_{k+1}$后，我们需要计算$u(x_{k+1})$和$u'(x_{k+1})$。为了控制精度，这里是通过泰勒展开式计算的。给定$\eqref{eq1}$式和在$x_k$处的$u(x_k)$和$u'(x_k)$，以及适当选取的整数$m \geq 2$,我们通过递推式计算出$u^{(3)}(x_k),u^{(4)}(x_{k}),\cdots,u^{(m)}(x_k)$，然后通过泰勒展开式估计$u(x_{k+1}),u'(x_{k+1})$,其中$h = x_{k+1} - x_{k}$。这样做的截断误差是$h^m$。
$$
u(x_{k+1}) = \sum_{k=0}^{m} \dfrac{u^{(k)}(x_k)}{k!} h^k + \epsilon \\
u’(x_{k+1}) = \sum_{k=1}^{m} \dfrac{u^{(k)} (x_k)}{(k-1)!} h^{k-1} + \tilde{\epsilon}
$$


根据论文所写，选择$m=30$可以达到16位精度，$m=60$可以达到34位精度，$m=120$可以达到80位精度。

### 简要叙述

#### 第一个根

给定一个初始值$x_s$,算法分两步计算大于或等于$x_s$的最小根$x_1$。

第一步，对于初始值$x_s$,通过Prufer变换$\eqref{eq2}$式计算出$\theta_{0} = \theta(x_s)$,求解初值$x(\theta_0) = x_s$和$\eqref{eq6}$式组成的微分方程：
$$
x(\theta_0) = x_s \\
\dfrac{dx}{d\theta} = -(\sqrt{\dfrac{r}{p}} + \dfrac{r'p - p'r + 2rq}{2rp} \dfrac{\sin(2\theta)}{2} )^{-1}
$$
使用Runge-Kutta法，在区间$\theta \in (\theta_0,-\dfrac{\pi}{2})$上计算，根据Prufer变换中的推导，在$\theta=-\dfrac{\pi}{2}$处的x值是多项式的一个根，得到$\tilde{x}_{1}$。

第二步，使用牛顿法，从$\tilde{x}_{1}$出发，计算得到高精度的根$x_{1}$。每一步都需要计算在$x_{i+1}$附近一点的$u$和$u'$，对于正交多项式，可以利用递推式$\eqref{eq8}$和$\eqref{eq9}$来计算函数值，代价与多项式的次数成正比。

#### 下一个根

给定$u$的一个根$x_i$,算法分两步计算得到$u$的下一个根$x_{i+1}$。

第一步，求解初值$x(\pi/2) = x_i$和$\eqref{eq6}$式组成的微分方程, 使用Runge-Kutta法在区间$\theta \in (\pi/2,-\pi/2)$上，在$\theta(x) = -\dfrac{\pi}{2}$处的x值是多项式的下一个根，得到$\tilde{x}_{i+1}$。

第二步，通过牛顿法提高$\tilde{x}_{i+1}$的精度。这里使用的$u$和$u'$值是通过泰勒展开估计的，$u$的高阶导数通过$\eqref{eq10}$来递推得到。

### 详细描述

算法1是用于给定第一个根，计算出多项式的所有根。

算法2和算法3用于计算多项式的第一个根。如果我们知道多项式的一个局部极值点，使用算法2，适用于勒让德和赫米特的偶数阶多项式。在泰勒展开收敛性不好的区域，使用算法3，通常只用于拉盖尔多项式。

#### 算法1

输入：多项式的常微分方程系数p,q,r;u的初始根$x_1$;大于等于$x_1$的根的数目N;在$x_1$处的导数$u'(x_1)$。

输出：N个根$x_1,x_2,\cdots,x_N$。

算法步骤：

令roots(1)=$x_1$,$ders(1)=u'(x_1)$。

do i=1,N=1

1. 给定初值$\theta_{0} = \pi/2$和$x(\theta_{0})=roots(i)$,在区间$[\pi/2,-\pi/2]$上用Runge-Kutta法求解$\eqref{eq6}$式。$x(-\pi/2)$即为roots(i+1)的初始估计值。
2. 使用递推式$\eqref{eq10}$计算得到u在roots(i)处的前m阶导数，从$u(x)=0,u'(x)=ders(i)$开始。
3. 使用牛顿法提高$roots(i+1)=x(-\pi/2)$的精度，其中$u$和$u'$是通过泰勒展开的前m项在roots(i)附近估计得到的。令ders(i+1) = u'(roots(i+1)）。

end do



#### 算法2

输入：多项式的常微分方程系数p,q,r;u的极值点$x_e$;在$x_e$处的函数值$u(x_e)$。

输出：满足$x_1 > x_s$的第一个最小根$x_1$。

算法步骤：

1. 给定初值$\theta_{0} = 0$和$x(\theta_{0})=x_e$,在区间$[0,-\pi/2]$上用Runge-Kutta法求解$\eqref{eq6}$式。$x(-\pi/2)$即为$x_1$的初始估计值。

2. 使用递推式$\eqref{eq10}$计算得到u在$x_e$处的前m阶导数，从$u(x_e),u'(x_e)=0$开始。

3. 使用牛顿法提高$x_{1}=x(-\pi/2)$的精度，其中$u$和$u'$是通过泰勒展开的前m项在**$x_e$附近**估计得到的。返回$x_1$。

#### 算法3

输入：函数表达式$u,u'$,多项式的常微分方程系数p,q,r;满足$x_{s} < x_{1}$的开始值$x_s$,$x_{1}$是u的最小根。

输出：u的最小根$x_1$。

算法步骤：

 1. 通过$\eqref{eq4}$式计算$\theta_{0} = \theta(x_s)$。

2. 给定初值$x(\theta_{0})=x_s$,在区间$[\theta_{0},-\pi/2]$上用Runge-Kutta法求解$\eqref{eq6}$式。$x(-\pi/2)$即为$x_1$的初始估计值。

3. 使用牛顿法提高$x_{1}=x(-\pi/2)$的精度，其中$u$和$u'$是通过给定的表达式计算得到。返回$x_1$。

### 测试

对函数$f(x) = x^{100}$，使用上述三种正交多项式得到的数值积分进行计算，根据n个点可以达到$2n-1$阶精度，采取51个点。计算中泰勒展开选择120项，牛顿迭代法的精度为80阶。

对于Lengendre多项式：

![image-20230508100052573](./assets/image-20230508100052573.png)

对于Hermite多项式:

![image-20230508100142722](./assets/image-20230508100142722.png)

对于Laguerre多项式：

![image-20230508100255297](./assets/image-20230508100255297.png)

Laguerre多项式积分的绝对误差很大。

先对$x^{60}$进行测试：

![image-20230508161645854](./assets/image-20230508161645854.png)

可以看到此时的误差已经无法达到80阶精度。

根据论文中的remark,对Laguerre多项式进行牛顿迭代法时，求高阶导数的公式为：
$$
\begin{aligned}
x u^{(k+2)} 
&= -(k + 1 -x) u^{(k+1)}  \\
&- (n - k) u^{(k)} 
\end{aligned}
$$
在x=0处是奇异点，导致使用泰勒展开项时，无法很好地逼近$x_{i+1}$,而是会更靠近$x_i$。

为了解决这一问题，我们改变计算流程，首先，使用Runge-Kutta法求出所有根的低阶精度逼近，不适用牛顿迭代法和泰勒展开。然后，与Lengendre和Hermite多项式不同，从最大根$x_{n}$开始，使用牛顿迭代法提高精度。

### 数值结果

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

