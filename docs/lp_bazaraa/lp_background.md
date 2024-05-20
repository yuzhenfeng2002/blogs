# 基础知识

这一章主要涉及线性代数、多面体的一些概念与常用定理，为后续学习线性规划的求解服务。

## 线性方程组

线性方程组可以用矩阵的形式表示为$Ax = b$，其中$A$为$m\times n$阶矩阵、$x$为$n$维未知数列向量、$b$为$m$维列向量。假设该方程组有解，则矩阵$A$的秩与增广矩阵$(A, b)$的秩满足$r(A) = r(A, b) = r \leq \min{\{ m, n \}}$。

我们可以通过初等行变换，将增广矩阵$(A, b)$变换为以下形式：
$$
(A, b) = \begin{bmatrix}
    A_1 & b_1\\
    A_2 & b_2
\end{bmatrix}
$$
其中，$A_1$为$r\times n$阶矩阵、且$A_1$的秩为$r$（因此，$A_1$的各行向量是线性无关的），$A_2$为$(m-r)\times n$阶矩阵，$b_1$为$r$维列向量，$b_2$为$(m - r)$维列向量。

由于$r(A) = r(A, b) = r$，且$r(A_1) = r$，因此满足$A_1x = b_1$的解$x$也满足$A_2x = b_2$，因此我们只需要考察$A_1x = b_1$即可。此时，$A_1$的各列向量中可以取出线性无关的$r$列（由于行秩等于列秩），组成新的矩阵，记作$B$，其余列组成矩阵$N$；$B$中各列对应的变量称为基变量、记作$x_B$，$N$中各列对应的变量称为非基变量、记作$x_N$。于是最初的线性方程组转化为了：
$$
(B, N)
\begin{bmatrix}
    x_B\\x_N
\end{bmatrix} = b_1
$$

也就是$Bx_B + Nx_N = b_1$。由于矩阵$B$为满秩矩阵，因此可逆，我们可以使用$x_N$来表示基向量，即$x_B = B^{-1}(b_1 - Nx_N)$。

从实际计算的过程来看，我们求解线性方程组应将增广矩阵通过高斯—若当消元法最终转化为下列形式（假设矩阵的前$r$列线性无关）：
$$
(A, b) \rightarrow \begin{bmatrix}
    I_r & Q & b'_1\\
    0 & 0 & b'_2\color{red}{(0)}
\end{bmatrix}
$$
此时，$x_B = b'_1 - Qx_N$。

## 多面体集合

**多面体集合**是有限个数的半空间的交集，常常表示为$\{x: Ax\leq b\}$，它是个特殊的凸集。线性规划的可行域便是一个多面体，我们通常考察$X = \{x: Ax \leq b \text{ and } x \geq 0\}$。

### 极点

在凸集中，有一些点不能够表示成两个不同点的严格凸组合，即：对于凸集$C$，如果$x = \lambda x_1 + (1 - \lambda)x_2\in C, \lambda \in (0, 1)$且$x_1, x_2\in C$，则$x = x_1 = x_2$，那么称$x$为凸集$C$的**极点**。

### 几何意义

对于**多面体集合**，假设我们考察的空间为$\mathbb{R}^n$，我们可以形象地把极点理解为线性无关的$n$个超平面的交点，也就是对于极点$\bar{x}$，有$n$个线性无关的约束是“紧”的（即只能取等号）。这个定义和之前的定义等价吗？

首先，我们说明若在点$x$处某约束是“紧”的，且$x$可以表示成$x_1, x_2$的严格凸组合，则$x_1, x_2$两点在该约束处也是“紧”的。根据条件，$Ax = A(\lambda x_1 + (1 - \lambda)x_2)\leq \lambda b + (1 - \lambda)Ax_2 \leq \lambda b + (1 - \lambda)b = b \leq b$，因此该式子中的众多不等号均取等，所以得到了我们的结论。

因此，若多面体集合中的一点$\bar{x}$，有$n$个线性无关的约束是“紧”的，则如果$\bar{x} = \lambda x_1 + (1 - \lambda)x_2\in X, \lambda \in (0, 1)$且$x_1, x_2\in X$，那么$\bar{x}, x_1, x_2$在这$n$个线性无关的约束上都是“紧”的。此时，有$Gx = g, x=\bar{x}, x_1, x_2$且$r(G) = n$，因此$\bar{x} = x_1 = x_2$。

反之，如果多面体集合中的一点$\bar{x}$处，有$r(r < n)$个线性无关的约束是“紧”的，此时有$G\bar{x} = g$，$G$为$(r\times n)$的矩阵、代表着“紧”的$r$个线性无关约束的系数。由于$r < n$，因此*一定*存在$d$，使得$Gd = 0$。且存在$\epsilon$，使得$x_1 = \bar{x} + \epsilon d\in X, x_2 = \bar{x} - \epsilon d\in X$，且$Gx_1 = g, Gx_2 = g$。由于$\bar{x} = 0.5x_1 + 0.5x_2$，所以$\bar{x}$不是极点。

形象地理解，$x$每被一个有效的约束“束缚”住，它的“自由度”就少了$1$，直到它的自由度变为$0$，也就成了极点。

### 极方向

给定一个集合$C$，它的**方向**指的是非零向量$d$，满足对于集合中的任意点$x$，集合$\{x + \lambda d: \lambda \geq 0\}\subseteq C$。对于多面体集合$X$，即：对于任意$x\in X$
$$
A(x + \lambda d) \leq b\\
x + \lambda d \geq 0
$$
由此可以得到：
$$
Ad\leq 0\\
d\geq 0\\
d \neq 0
$$

将方向$d$标准化（$1^Td = 1$），得到集合$D = \{d: Ad\leq 0\text{ and } d\geq 0\text{ and } 1^Td = 1\}$，我们称为多面体集合$X$的退缩方向集。

对于一个集合$C$所有方向$D$，如果$d = \lambda d_1 + (1 - \lambda)d_2\in D, \lambda > 0$且$d_1, d_2\in D$，则$d = d_1 = d_2$，那么称$d$为集合$C$的**极方向**。集合$D$的极点对应着集合$X$的极方向。
<!-- 对于凸锥（满足$\forall x \in C \text{ and } \lambda \geq 0, \lambda x \in C$的凸集）来说，知道了它的极方向就知道了它。 -->

如下图所示，解释了多面体的极点、极方向、退缩锥以及它们之间的关系。

![](docs/lp_bazaraa/lp_poly_extreme.png)

## 表示定理

令$X = \{x: Ax \leq b \text{ and } x\geq 0\}$为一非空的多面体集合，则存在极点$x_1, \cdots, x_k$，其中$1\leq k < \infty$；集合$X$存在极方向当且仅当$X$无界，若$X$无界，则存在极方向$d_1, \cdots, d_l$，其中$1\leq l < \infty$；$\bar{x}\in X$当且仅当$\bar{x}$可以表示为$x_1, \cdots, x_k$的凸组合与$d_1, \cdots, d_l$的非负线性组合之和。

### 证明

若证明了$1\leq k < \infty$，那么对于集合$D = \{d: Ad\leq 0\text{ and } d\geq 0\text{ and } 1^Td = 1\}$，其极点的个数$l$也满足$1\leq l < \infty$，对应的极方向的个数也就$1\leq l < \infty$。因此，对于极点和极方向的存在性与个数，我们只需要证明极点的个数满足$1\leq k < \infty$即可。

在多面体集合中，选择一点$\bar{x}\in X$，记在该点处的线性无关的约束个数为$r$。若$r = n$，则该点为极点，$k\geq 1$；若$r < n$，则根据定义有$\bar{x} = \lambda y_1 + (1 - \lambda)y_2, \lambda\in (0, 1), y_1, y_2\in X, y_1\neq y_2$。此时，令$d = y_2 - y_1$，则有${y}_{1}=\bar{{x}}-(1-\lambda) {d}$且${y}_{2}=\bar{{x}}+\lambda {d}$，这说明从$\bar{x}$沿着$d$与$-d$方向均*可*进行“游走”。令$\bar{\gamma} = \max\{\gamma: \bar{x} - \gamma d\in X\}$，因为$X\subseteq \{x: x\geq 0\}$，所以$\bar{\gamma} < \infty$。令$\bar{y_1} = \bar{x} - \bar{\gamma} d$。可以直观地想象，$\bar{y_1}$一定是碰到了“约束”才不能够继续沿着$-d$方向“游走”。所以，在$\bar{y_1}$点处的线性无关的约束个数至少为$r + 1$。令$\bar{x} = \bar{y_1}$，重复该步骤，直到在$\bar{x}$处的线性无关的约束个数为$n$，即为极点。综上所述，$k \geq 1$。又由于从$m + n$（连同$x\geq 0$这$n$条约束）条约束中选择$n$条约束只有有限种选择方式，因此$k < \infty$。

若$\bar{x}$可以表示为$x_1, \cdots, x_k$的凸组合与$d_1, \cdots, d_l$的非负线性组合之和，很容易验证$A\bar{x}\leq b\text{ and } x\geq 0$，即$\bar{x}\in X$。若$\bar{x}\in X$，它如何表示成$x_1, \cdots, x_k$的凸组合与$d_1, \cdots, d_l$的非负线性组合之和的形式呢？

我们考察集合
$$
\bar{X} = X\cap\{x:1^Tx\leq M\}
$$
里的一点$\bar{x}$，其中$M$足够大，此时极点集合$\bar{S_p} = \{x_1, \cdots, x_k, x_{k + 1}, \cdots, x_{k+u}\}$，我们首先证明$\bar{x}$可以用集合$\bar{X}$的极点的凸组合表示。

如果$\bar{x}\in \bar{S_p}$，则结论成立；如果$\bar{x}\not \in \bar{S_p}$，设$Gx = g$为$\bar{x}$满足的“紧”的约束，由于$r(G)\leq n - 1$，所以一定存在$d\neq 0, Gd = 0$。与上面类似，将$\bar{x}$沿着$d$“游走”，直到游走至极点$\bar{y_1}$。

再令$\bar{y_2} = \bar{x} + \bar{\gamma_2}(\bar{x} - \bar{y_1})$，其中$\bar{\gamma_2} = \max\{\gamma: \bar{x} + \gamma(\bar{x} - \bar{y_1})\in \bar{X}\}$。此时，$\bar{x} = \frac{1}{1 + \bar{\gamma_2}}(\gamma_2 \bar{y_1} + \bar{y_2})$。注意到$G\bar{y_2} = g$，同时在$\bar{y_2}$处的线性无关的约束个数至少比$\bar{x}$增加了$1$。因此，若$\bar{y_2}$是$\bar{X}$的极点，则结论得证；若不是，则继续将$\bar{y_2}$表示成两个点的凸组合，与以上步骤类似。最终，可以得到：
$$
\bar{x} = \sum\limits_{j = 1}^{k + u}\delta_jx_j, \sum\limits_{j = 1}^{k + u}\delta_j = 1\text{ and } \delta\geq 0
$$

如果$\delta_j = 0, j = k + 1, \cdots, k + u$，则定理得证。否则，由于$x_j, j = k + 1, \cdots, k + u$是约束$1^Tx\leq M$和其他$n - 1$条线性无关的约束形成的极点。从中取极点$x_v$，令$G$表示该点满足的$n - 1$条约束，则$Gd = 0$有非零解$d$。$x_v$可沿着$d$游走至它的邻接点$x_{i(v)}$，则$\bar{d} = \frac{x_v - x_{i(v)}}{1^T(x_v - x_{i(v)})}$是集合$X$的方向（因为没有了$1^Tx\leq M$这条约束）。对于多面体集合$D$，$\bar{d}$满足了其中的$n - 1$条约束，又满足了等式约束$1^Td = 1$，因此$\bar{d}$是集合$X$的极方向，所以$x_v = x_{i(v)} + \theta_v\bar{d}, \theta_v = 1^T(x_v - x_{i(v)})$。

最终，我们得到了：
$$
\bar{{x}}=\sum_{j=1}^{k} \delta_{j} {x}_{j}+\sum_{v=k+1}^{k+u} \delta_{v} {x}_{i(v)}+\sum_{v=k+1}^{k+u} \delta_{v} \theta_{v} {d}_{j(v)}
$$
即有
$$
\bar{{x}}=\sum_{j=1}^{k} \lambda_{j} {x}_{j} +\sum_{j=1}^{l} \mu_{j} {d}_{j} \\
\sum_{j=1}^{k} \lambda_{j}=1 \\
\lambda_{j}  \geq 0, \quad j=1, \ldots, k \\
\mu_{j} \geq 0, \quad j=1, \ldots, l
$$

注意，我们只有不超过$\min \{n + 1, k + l\}$个$\lambda, \mu$是非零的。原因在于，上式中若以$(\lambda, \mu)$为变量，则等式的秩$r\leq \min \{n + 1, k + l\}$；同时，上述约束也构成了一多面体，也至少有一个极点，在该极点处有至少$k + l$个线性无关的约束是“紧”的。因此，还需有$k + l - r$个紧约束，即这些变量为$0$。