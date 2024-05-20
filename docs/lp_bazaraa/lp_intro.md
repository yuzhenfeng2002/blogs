# 线性规划介绍

## 线性规划问题

我们称形如
$$
\begin{array}{lrrrrrrr}
\text {minimize} & c_{1} x_{1} & + & c_{2} x_{2} & +\cdots+ & c_{n} x_{n} & & \\
\text {subject to} & a_{11} x_{1} & + & a_{12} x_{2} & +\cdots+ & a_{1 n} x_{n} & \geq & b_{1} \\
& a_{21} x_{1} & + & a_{22} x_{2} & +\cdots+ & a_{2 n} x_{n} & \geq & b_{2} \\
& \vdots & & \vdots & +\cdots+ & \vdots & & \vdots \\
& a_{m 1} x_{1} & + & a_{m 2} x_{2} & +\cdots+ & a_{m n} x_{n} & \geq & b_{m} \\
& x_{1}& , & x_{2}& ,\cdots, & x_{n} & \geq & 0
\end{array}
$$
的问题为**线性规划**（Linear Programming）问题。其中，$c_{1} x_{1} +  c_{2} x_{2} + \cdots + c_{n} x_{n}$称为**目标函数**，这是我们需要优化（这里是最小化）的目标；$x_j, j = 1, \cdots, n$称为**决策变量**，这是我们需要决策的内容；$c_j, j = 1, \cdots, n$称为成本系数，这代表着每增加单位数量的$x_j$所需要花费的成本；不等式$\sum_{j = 1}^n a_{ij}x_j \geq b_i, i = 1, \cdots, m$称为**约束**；由系数$\mathbf{A} = (a_{ij})_{m\times n}$组成的矩阵称为约束矩阵。满足所有约束不等式的解$(x_1, \cdots, x_n)$称为一个**可行解**，所有可行解的集合构成了**可行域**，求解该线性规划问题即在可行域中寻找可行解使得目标函数最小。

很多其他形式的线性规划问题可以与上述形式相互转化，包括：
- 等式与不等式约束之间相互转化
  - 等式约束$\sum_{j = 1}^n a_{ij}x_j = b_i$等价于$\sum_{j = 1}^n a_{ij}x_j \geq b_i\text{ and }\sum_{j = 1}^n a_{ij}x_j \leq b_i$；
  - 不等式约束$\sum_{j = 1}^n a_{ij}x_j \geq b_i$等价于$\sum_{j = 1}^n a_{ij}x_j - x_{n + i} = b_i$，同理$\sum_{j = 1}^n a_{ij}x_j \leq b_i$等价于$\sum_{j = 1}^n a_{ij}x_j + x_{n + i} = b_i$；
- 决策变量的非负性
  - 无约束变量$x_j$等价于两个非负变量之差，即$x_j = x'_j - x''_j, x'_j, x''_j\geq 0$；
  - 若有多个无约束变量$x_1, \cdots, x_k$，则令$x_j = x'_j - x'', x''\geq 0$即可（只需要一个$x''$变量）；
- 最小化目标函数
  - 若线性规划的目标为最大化目标函数，则可以通过$\max \sum_{j=1}^n c_jx_j = -\min -\sum_{j=1}^n c_jx_j$进行转换。

线性规划的几种形式较为常见，包括：
- 标准型（常见于单纯形法）：所有约束均为等式，且所有变量均非负；
- 权威型（常见于对偶性）：
  - 对于最小化问题，所有约束均为$\geq$形式，且所有变量均非负；
  - 对于最大化问题，所有约束均为$\leq$形式，且所有变量均非负。

我们可以使用矩阵形式更加方便地表达线性规划问题，例如对于：
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \sum\limits_{j=1}^n c_jx_j \\
\text {subject to} & \sum\limits_{j = 1}^n a_{ij}x_j & = & b_i,&  i=1, \cdots, m \\
& x_j & \geq & 0,&  j=1, \cdots, n
\end{array}
$$
令
$$
\mathbf{c}=\left[\begin{array}{c}
c_{1} \\ c_{2} \\\vdots \\ c_{n}
\end{array}\right], 
\mathbf{x}=\left[\begin{array}{c}
x_{1} \\ x_{2} \\\vdots \\ x_{n}
\end{array}\right], \mathbf{b}=\left[\begin{array}{c}
b_{1} \\ b_{2} \\\vdots \\ b_{m}
\end{array}\right], \mathbf{A}=\left[\begin{array}{cccc}
a_{11} & a_{12} & \ldots & a_{1 n} \\
a_{21} & a_{22} & \ldots & a_{2 n} \\
\vdots & \vdots & & \vdots \\
a_{m 1} & a_{m 2} & \ldots & a_{m n}
\end{array}\right]
$$
于是，该线性规划问题可以描述为：
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \mathbf{c}^T\mathbf{x} \\
\text {subject to} & \mathbf{Ax} & = & \mathbf{b}\\
& \mathbf{x} & \geq & \mathbf{0}
\end{array}
$$
若将矩阵$\mathbf{A}$的每列作为一个向量，即记$\mathbf{A} = \left[\mathbf{a}_{1}, \mathbf{a}_{2}, \cdots, \mathbf{a}_{n}\right]$，则原问题可以描述为：
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \sum\limits_{j=1}^n c_jx_j \\
\text {subject to} & \sum\limits_{j = 1}^n \mathbf{a}_{j}x_j & = & \mathbf{b}\\
& x_j & \geq & 0,&  j=1, \cdots, n
\end{array}
$$

## 几何法求解

对于维度较低（二维或三维）的线性规划问题，可以使用几何法进行求解。

$\sum\limits_{j = 1}^na_{ij}x_j  =  b_i$表示的是一个（超）平面（二维则是直线），而$\sum\limits_{j = 1}^na_{ij}x_j  \geq  b_i$则是（超）平面朝向$(a_{i1}, \cdots, a_{in})$的区域。通过约束，可以画出可行域。对于目标函数，令$\sum\limits_{j=1}^n c_jx_j = z$，这也是一个（超）平面。那么，若原问题是最小化问题，则将法向量为$\mathbf{c}$的平面沿着逆法向量（即$-\mathbf{x}$）移动至可行域的边界处，此处（若存在）即为最优解；若原问题是最大化问题，则沿法向量方向移动。如下图所示。

<img src="docs/lp_bazaraa/lp_geo_sol.png" alt="几何法求解" width="200"/>

从几何法可以看到，线性规划问题有以下几种结果：
- 有唯一最优解。无论是可行域有界或无界，都有可能出现这种情况；
  - <img src="docs/lp_bazaraa/lp_geo_unique.png" alt="唯一最优解" width="400"/>
- 有无穷多最优解，但最优解存在（$\mathbf{c}^T\mathbf{x}^* < \infty$）。同样，可行域有界或无界，都有可能出现这种情况；
  - <img src="docs/lp_bazaraa/lp_geo_alter.png" alt="无穷多最优解" width="400"/>
- 最优解无界（$\mathbf{c}^T\mathbf{x}^* = \infty$）；
  - <img src="docs/lp_bazaraa/lp_geo_unbounded.png" alt="最优解无界" width="200"/>
- 可行域为空。

我们也可以基于另一视角考察线性规划的几何意义。对于问题：
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \sum\limits_{j=1}^n c_jx_j \\
\text {subject to} & \sum\limits_{j = 1}^n \mathbf{a}_{j}x_j & = & \mathbf{b}\\
& x_j & \geq & 0,&  j=1, \cdots, n
\end{array}
$$
若向量$\mathbf{b}$在由$\mathbf{a}_{1}, \mathbf{a}_{2}, \cdots, \mathbf{a}_{n}$围成的“锥”中，则该问题可行（可行域非空）。同时，令$z = \sum\limits_{j=1}^n c_jx_j$，求解该线性规划问题即在向量
$$\left[\begin{array}{l}c_{1} \\ \mathbf{a}_{1}\end{array}\right],\left[\begin{array}{l}c_{2} \\ \mathbf{a}_{2}\end{array}\right], \cdots,\left[\begin{array}{l}c_{n} \\ \mathbf{a}_{n}\end{array}\right]$$
张成的“锥”中寻找使得$\left[\begin{array}{l}z \\ \mathbf{b}\end{array}\right]$中$z$最小的情况。

例如，对于线性规划问题：
$$
\begin{array}{lrrrrr}
\operatorname{minimize} & -2 x_{1} & - & 3 x_{2} \\
\text {subject to} & x_{1} & + & 2 x_{2} & \leq & 2\\
& x_1 & , & x_2 & \geq & 0
\end{array}
$$
将不等式约束转化为等式约束，得到：
$$
\left[\begin{array}{r}
-2 \\ 1
\end{array}\right] x_{1}+\left[\begin{array}{r}
-3 \\ 2
\end{array}\right] x_{2}+\left[\begin{array}{l}
0 \\ 1
\end{array}\right] x_{3}=\left[\begin{array}{l}
z \\ 2
\end{array}\right]
$$
其几何意义如下图所示。

<img src="docs/lp_bazaraa/lp_req_space.png" alt="几何法求解" width="300"/>

其中，阴影部分即为向量围成的“锥”，虚线表示的是向量$\mathbf{b}$，在阴影部分里的虚线段上，寻找到的使得$z$最小的点即为最优解。