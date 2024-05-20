# 单纯形法

## 极点与最优性

通过第一章的几何法，我们或许可以注意到：如果存在最优解，那么一定存在极点，在那里能够取到最优解。我们上一章中学习了表示定理，可以通过表示定理证明该结论。

对于线性规划问题
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \mathbf{c}^T\mathbf{x} \\
\text {subject to} & \mathbf{Ax} & = & \mathbf{b}\\
& \mathbf{x} & \geq & \mathbf{0}
\end{array}
$$
而言，其可行域中的任意一点$\bar{x}$，其满足：
$$
\bar{{x}}=\sum_{j=1}^{k} \lambda_{j} {x}_{j} +\sum_{j=1}^{l} \mu_{j} {d}_{j} \\
\sum_{j=1}^{k} \lambda_{j}=1 \\
\lambda_{j}  \geq 0, \quad j=1, \ldots, k \\
\mu_{j} \geq 0, \quad j=1, \ldots, l
$$
其中，可行域的极点为$x_1, \cdots, x_k$，极方向为$d_1, \cdots, d_l$。则该线性规划问题可以等价转化为
$$
\min \sum_{j=1}^{k} \lambda_{j} (c^T{x}_{j}) +\sum_{j=1}^{l}\mu_{j} (c^T{d}_{j}) \\
\text{s.t. }
\begin{cases}
    \sum_{j=1}^{k} \lambda_{j}=1 \\
    \lambda_{j}  \geq 0, \quad j=1, \ldots, k \\
    \mu_{j} \geq 0, \quad j=1, \ldots, l
\end{cases}
$$

注意到$\mu_j$可以无限大，因此解存在下列几种情形：
- 如果存在$(c^T{d}_{j})< 0$，则该问题的最优解无下界（即不存在最优解）；
- 否则，取$\mu_j$满足$\mu_j(c^T{d}_{j}) = 0$。此时，记$j^\star = \arg\min(c^T{x}_{j})$，则令$\lambda_{j^\star} = 1$。此时
  $$
  \bar{{x}}=\sum_{j=1}^{k} \lambda_{j} {x}_{j} +\sum_{j=1}^{l} \mu_{j} {d}_{j} = x_{j^\star} + \sum_{c^T{d}_{j} = 0} \mu_{j} {d}_{j}
  $$
  取最优解。

显然，当取$\mu_j = 0$时，$\bar{x}$即为极点。另外，若$c^T{x}_{j}$取最小时$j$不唯一，则最优解可以表示为
$$
\bar{{x}}=\sum_{c^T{x}_{j} = c^T{x}_{j^\star}} \lambda_{j} {x}_{j} + \sum_{c^T{d}_{j} = 0} \mu_{j} {d}_{j}\\
\sum_{c^T{x}_{j} = c^T{x}_{j^\star}} \lambda_{j} = 1\\
\mu_{j}\geq 0
$$

## 基可行解

极点是几何视角下的词，极点是和代数视角下的基可行解对应的（但并非一一对应）。下面给出基可行解的定义。

对于约束$Ax = b, x\geq 0$（其中$A$是$m\times n$的矩阵，$b$为$m$维向量），假设$r(A, b) = r(A) = m$。此时，$A$的各列向量中可以取出线性无关的$m$列（由于行秩等于列秩），组成新的矩阵，记作$B$，其余列组成矩阵$N$；$B$中各列对应的变量称为基变量、记作$x_B$，$N$中各列对应的变量称为非基变量、记作$x_N$。于是最初的线性方程组转化为了：
$$
(B, N)
\begin{bmatrix}
    x_B\\x_N
\end{bmatrix} = b_1
$$
也就是$Bx_B + Nx_N = b$。由于矩阵$B$为满秩矩阵，因此可逆，我们可以使用$x_N$来表示$x_B$，即$x_B = B^{-1}(b - Nx_N)$。

令$x_N = 0$，则解$x = (x_B, x_N)^T$
$$
x_B = B^{-1}b\\
x_N = 0
$$
称为**基解**；若$x_B\geq 0$，则$x$称为**基可行解**（因为它可行了）；如果$x_B > 0$，则$x$称为**非退化**基可行解；反之，若$x_B$中有值为$0$，则$x$称为**退化**基可行解。几何视角下的极点是和代数视角下的基可行解对应的，上一章我们把定义极点为线性无关的$n$个超平面的交点，若在这一点有超过$n$个线性无关的约束是“紧”的，则称这个极点是退化极点，对应的基可行解也是退化基可行解。

为什么说极点和基可行解是对应的呢？一方面，一个点是极点意味着它对于$n$个线性无关的约束是“紧”的，但是$m \leq n$，因此还需要$p = n - m$个“紧”的约束，因此需要$x_N = 0$作为剩下的$p$个约束，而极点$x = (x_B, x_N)^T$又肯定是可行的。另一方面，基可行解满足$x_B = B^{-1}b, x_N = 0$，满足$Ax = b, x_N = 0$，这是$n$个线性无关的“紧”约束，因此符合极点的定义。

前面我们说极点和基可行解只是对应、而不是一一对应的关系，问题在于一个极点可能对应多个基可行解，这种情况发生在退化的情况。此时，在该极点处有多于$n$个的线性无关的“紧”约束，那么$x_B$中就有值也为$0$（本不该为$0$，因为约束$x_B\geq 0$本不应该是“紧”的），于是我们只要把$x_B$中不为$0$的部分对应的列作为$B$的一部分、其他部分任意选择，即可得到对应该极点的基可行解。

但是，并不是所有的退化极点都对应多个基可行解的。我们讨论下面这种情况，$Ax = b$“隐含”着$x_k = 0$这一条约束，那么一方面若$B$中包括$x_k$对应的这一列，则解$x = (x_B, x_N)^T$满足$Ax = b, x_N = 0, x_k = 0$这些约束，是退化解；另一方面，由于$Ax = b$“隐含”着$x_k = 0$，那么矩阵$(A, b)$可以通过初等行变换转化为形如
$$
\begin{matrix}
    \text{The 1st row: } & 1 & \cdots &  & & \cdots & \bar{b_1}\\
    \vdots & \vdots & \ddots\\
    \text{The kth row: } & 0 & \cdots & 1 & 0 & \cdots & 0\\
    \text{The (k+1)th row: } & 0 & \cdots & 0 & 1 & \cdots & \bar{b_{k + 1}}\\
    \vdots
\end{matrix}
$$
可以看到，若排除$x_k$对应的那一列，矩阵的秩只有$m - 1$，因此$x_k$只能作为基变量，基可行解不可以有多个。

从上面的论述中，我们知道了如果一个线性规划问题有有界的最优解，那么它一定存在能取到最优解的极点，而极点又是和基可行解对应的，一个基可行解对应一个极点。那么，很自然的想法是对于所有的基可行解，算出目标函数值进行比较，即可求得最优的那个基可行解。但这样的想法有三个缺点：
- 基解数量很多，最多有$\left(\begin{matrix}n\\m\end{matrix}\right)$个基可行解；
- 无法判断最优解是否无界；
- 对于不可行的情况，也需要算出所有的基解来进行判断。

因此，单纯形法诞生了。单纯形法的使命就是“聪明”地进行遍历。

## 单纯形法

### 最优性判别

人的一生，或许就像在一片汪洋中寻找着埋藏在最深处的宝藏，但是站在这片汪洋的某一点，往往不知道哪里才是最深处，也不知道是否已经到达了最深处。而这对于线性规划来说，很容易，只可惜人生的约束不是线性的。

对于线性规划问题
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \mathbf{c}^T\mathbf{x} \\
\text {subject to} & \mathbf{Ax} & = & \mathbf{b}\\
& \mathbf{x} & \geq & \mathbf{0}
\end{array}
$$
假设$r(A, b) = r(A) = m$，其中$A$是$m\times n$的矩阵，$b$为$m$维向量。对于一组基，我们有$Ax = Bx_B + Nx_N = b$，因此
$$
x_B = B^{-1}b - B^{-1}Nx_N
$$
此时，目标函数值
$$
z = c^Tx = c^T(x_B, x_N)^T\\
= c_B^Tx_B + c_N^Tx_N \\
= c_B^T(B^{-1}b - B^{-1}Nx_N) + c_N^Tx_N \\
= c_B^TB^{-1}b + (c_N^T - c_B^TB^{-1}N)x_N
$$
考虑此时的基可行解$x = (x_B, x_N)^T = (B^{-1}b, 0)^T$，在该解处的目标值为$z_0 = c_B^TB^{-1}b$，那么原线性规划问题可以等价转化为：
$$
\begin{array}{lrrrr}
\operatorname{minimize} & (\mathbf{c}_N^T - \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{N})\mathbf{x}_N & + & z_0 \\
\text {subject to} & \mathbf{B}^{-1}\mathbf{N}\mathbf{x}_N & + & \mathbf{x}_B & = & \mathbf{B}^{-1}\mathbf{b}\\
& \mathbf{x}_N & , & \mathbf{x}_B & \geq & \mathbf{0}
\end{array}
$$
注意上述过程中只涉及了变量替换。上述形式可以继续等价转化为：
$$
\begin{array}{lrrrr}
\operatorname{minimize} & (\mathbf{c}_N^T - \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{N})\mathbf{x}_N & + & z_0 \\
\text {subject to} & \mathbf{B}^{-1}\mathbf{N}\mathbf{x}_N & \leq & \mathbf{B}^{-1}\mathbf{b}\\
& \mathbf{x}_N & \geq & \mathbf{0}
\end{array}
$$
从上式可以看到，当$\mathbf{c}_N^T - \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{N}\geq 0$时，当前基可行解就是最优的（不管$x_N$取什么值，目标函数值都不会比$x_N = 0$时小），这就像你在茫茫大海中知道这里就是世界最深的地方。

### 几何意义

下面，我们先以$p = n - m = 2$为例，探究一下单纯形法迭代的几何意义。我们考察如下图所示的情形：

<img src="docs/lp_bazaraa/lp_simplex_iter.png" alt="几何意义" width="300"/>

其中，$(x_3, x_4, x_5, x_6)$为基变量，$(x_1, x_2)$为非基变量，各条约束对应着图中的以直线（直线对应着基变量等于$0$）为界的半空间。当非基变量为$0$时，即$x_1 = x_2 = 0$，对应着图中原点$v_1$。

由于$-\bar{c} = \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{N} - \mathbf{c}_N^T > 0$，因此不管是增大$x_1$还是$x_2$，都可以带来目标函数值的减小。此时，我们不妨增大$x_2$（这里不妨随意选择，但在单纯形法中，我们不是随意选的），每增加单位$x_2$，我们获得的目标函数的减少量为$\bar{c}_2$。

但是，由于约束的存在，$x_2$不能无限增加，它会一直增加到$x_3 = 0$，再增加下去$x_3$便要为负了，这是不允许的。于是我们到了$v_2$点，此时$x_1 = x_3 = 0$，这两个变量为基变量。

接着我们看到，增加$x_1$的方向和$-\bar{c}$的夹角仍为锐角，我们可以增加$x_1$以减少目标函数值，一直到$x_4$或$x_5$为$0$，即到达了$v_3$点。$v_3$点是一个退化极点，我们不妨选择$x_3$和$x_4$为非基变量，但是此时即使$x_3$只增加一点点，$x_5$就会变成负的，所以我们再将$x_3$变为基变量，$x_5$变为非基变量，再增加$x_4$至$v_4$点，$x_5$和$x_6$成为非基变量。在这里，无论是增加$x_5$还是增加$x_6$，都是和$-\bar{c}$呈钝角，没有办法使得目标函数值减少了，就到达了最优的极点处。

因此，单纯形法的每次迭代，就是视当前的基对应的极点为原点，在原点处看可不可以向各个轴的正方向进行增加；如果有，则选择一根轴，一直增加到到达另一个极点处的时候，再把这个极点视为原点，开始下一轮迭代。而这里的“单纯形”指的是在$\mathbb{R}^P$空间中，原点和各个轴上与它相邻的极点这$p + 1$个*不在同一个超平面上点*组成的多面体。

### 代数描述

下面，我们把上述过程描述成代数形式。我们已经把一个线性规划问题转化成了：
$$
\begin{array}{lrrrr}
\operatorname{minimize} & (\mathbf{c}_N^T - \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{N})\mathbf{x}_N & + & z_0 \\
\text {subject to} & \mathbf{B}^{-1}\mathbf{N}\mathbf{x}_N & + & \mathbf{x}_B & = & \mathbf{B}^{-1}\mathbf{b}\\
& \mathbf{x}_N & , & \mathbf{x}_B & \geq & \mathbf{0}
\end{array}
$$

我们将选择一个能使目标函数变小的非基变量，让它增加。直觉上看，我们选择$x_k$，其中$k = \arg\min\limits_i \{(\mathbf{c}_N^T - \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{N})_i: (\mathbf{c}_N^T - \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{N})_i < 0\}$。但是，$x_k$一般不能无限增大（最优解无界的情况除外），因为$x_B$必须是非负的。我们假设这个非基变量从$0$变为$x_k, x_k > 0$，则此时

目标函数值为
$$
z = z_0 - \left[\mathbf{c}_B^T(\mathbf{B}^{-1}\mathbf{N})_k - (\mathbf{c}_N)_k\right]x_k
$$

$x_B$变为了
$$
\mathbf{x}_B = \mathbf{B}^{-1}\mathbf{b} - (\mathbf{B}^{-1}\mathbf{N})_k x_k
$$

记
$$
\mathbf{B}^{-1}\mathbf{N} = (\mathbf{\bar{a}}_{N_1}, \cdots, \mathbf{\bar{a}}_{N_p}) = \begin{pmatrix}
    \mathbf{y}_{1}\\ \vdots \\ \mathbf{y}_{m}
\end{pmatrix}
$$
于是
$$
z = z_0 - \left(\mathbf{c}_B^T\mathbf{\bar{a}}_{k} - {c}_k\right)x_k\\
\mathbf{x}_B = \mathbf{\bar{b}} - \mathbf{\bar{a}}_{k} x_k
$$

由于$\mathbf{x}_{B}\geq 0$，因此
$$
\mathbf{x}_{B} = \begin{pmatrix}
    \bar{b_1}\\ \vdots \\ \bar{b_m}
\end{pmatrix} - 
\begin{pmatrix}
    {y_{1k}}\\ \vdots \\ {y_{mk}}
\end{pmatrix}x_k \geq 0
$$

取$r = \arg\min\limits_{i}\{\frac{\bar{b_i}}{y_{ik}}: y_{ik} > 0\}$，则$x_k$最多可增加至$x_k = \frac{\bar{b_r}}{y_{rk}}$，此时，
$$
\mathbf{x}_{B} = \begin{pmatrix}
    \bar{b_1}\\ \vdots \\ \bar{b_m}
\end{pmatrix} - 
\begin{pmatrix}
    {y_{1k}}\\ \vdots \\ {y_{mk}}
\end{pmatrix}\frac{\bar{b_r}}{y_{rk}} = \begin{pmatrix}
    \bar{b_1} - \frac{y_{1k}}{y_{rk}}\bar{b_r}\\
    \vdots\\
    0\\
    \vdots\\
    \bar{b_m} - \frac{y_{mk}}{y_{rk}}\bar{b_r}
\end{pmatrix}\\
\mathbf{x}_{N} = \begin{pmatrix}
    0\\
    \vdots\\
    \frac{\bar{b_r}}{y_{rk}}\\
    \vdots\\
    0
\end{pmatrix}\\
z = z_0 - \left(\mathbf{c}_B^T\mathbf{\bar{a}}_{k} - {c}_k\right)\frac{\bar{b_r}}{y_{rk}}
$$

若不考虑退化的情况（即$\bar{b_r} = 0$），则可以看到第$r$个基变量变为$0$，而第$k$个基变量从$0$变为$\frac{\bar{b_r}}{y_{rk}}$，实现了换基迭代操作。同时，目标函数值减少了$\left(\mathbf{c}_B^T\mathbf{\bar{a}}_{k} - {c}_k\right)\frac{\bar{b_r}}{y_{rk}}$。由于我们只有有限的基，且目标函数值一直在减少，因此这样的迭代过程不可能遇到重复的基（暂时不考虑退化的情况），所以迭代一定能够终止。

这样的结果也有一定经济意义。我们发现，$x_k$每增加$1$单位，基变量$x_{Bi}$就要减少$y_{ik}$，而$z = z_0 - \left(\mathbf{c}_B^T\mathbf{\bar{a}}_{k} - {c}_k\right)x_k = z_0 - \sum\limits_{i=1}^{m} (c_{Bi}y_{ik})x_k + {c}_kx_k$，因此我们可以视$\sum\limits_{i=1}^{m} (c_{Bi}y_{ik})$为增加$x_k$可以节省的基变量的成本，而$c_k$是增加$x_k$自身的成本。所以，也只有$\sum\limits_{i=1}^{m} (c_{Bi}y_{ik}) - x_k > 0$时增加$x_k$才会降低总成本。

### 迭代终止

通过以上的描述，我们知道单纯形法的每一次迭代都有进基和出基的操作，进基和出基的变量选择需要满足各自的条件，其中：
- 进基：$x_k$满足$\mathbf{c}_B^T\mathbf{\bar{a}}_{k} - {c}_k>0$；
- 出基：$x_r$满足$r = \arg\min\limits_{i}\{\frac{\bar{b_i}}{y_{ik}}: y_{ik} > 0\}$。

既然如此，就出现两种迭代终止的情况。

第一，是对于所有的非基变量$x_j$，$\mathbf{c}_B^T\mathbf{\bar{a}}_{j} - {c}_j\leq 0$，此时我们得到了最优解。这时，如果对于所有的非基变量$x_j$，$\mathbf{c}_B^T\mathbf{\bar{a}}_{j} - {c}_j < 0$，则最优解是唯一的；否则，若存在$\mathbf{c}_B^T\mathbf{\bar{a}}_{k} - {c}_k = 0$，则增加$x_k$不会带来最优解的减少，直到$x_k$不能再增大则得到了另一个基解。

第二，是对于所有的$i$，$y_{ik}\leq 0$，这意味着$x_k$能无限增大，最优解也就无界。事实上，取
$$
(\mathbf{x}_B,\mathbf{x}_N)^T = \begin{pmatrix}
    \bar{b_1}\\ \vdots \\ \bar{b_m}\\0\\ \vdots\\0\\ \vdots\\0
\end{pmatrix} - 
\begin{pmatrix}
    {y_{1k}}\\ \vdots \\ {y_{mk}}\\0\\ \vdots\\1\\ \vdots\\0
\end{pmatrix}x_k, x_k \geq 0
$$
即可，此时目标函数值为
$$
\mathbf{cx} = \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{b} - (\mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{{a}}_{k} - {c}_k)x_k
$$

### 单纯形法

对于线性规划问题
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \mathbf{c}^T\mathbf{x} \\
\text {subject to} & \mathbf{Ax} & = & \mathbf{b}\\
& \mathbf{x} & \geq & \mathbf{0}
\end{array}
$$
其中，$\mathbf{A}$是$m\times n$的矩阵，$b$为$m$维向量，且$r(\mathbf{A}) = m$，单纯形法的步骤如下：

1. 初始化。寻找基$\mathbf{B}$；
2. 求$\mathbf{B}^{-1}$，进而求$\mathbf{x}_B = \mathbf{B}^{-1}\mathbf{b}$，$\sigma_j = \mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{{a}}_{j} - {c}_j$，$\mathbf{\bar{A}} = \mathbf{B}^{-1}\mathbf{A}$；
3. 若$\sigma_j \leq 0$，则得到最优解，终止；否则，取$k = \arg\max\limits_j \{\mathbf{c}_B^T\mathbf{B}^{-1}\mathbf{{a}}_{j} - {c}_j\}$（这种取法称为Dantzig规则），以$x_k$为进基变量；
4. 若$\mathbf{a}_k\leq 0$，则最优解无界，终止；否则，取$r = \arg\min\limits_{i}\{\frac{\bar{b_i}}{y_{ik}}: y_{ik} > 0\}$，以$x_{Br}$为出基变量；
5. 换基，继续第2步。

我们可以用表格的形式来描述单纯形法迭代的过程，我们将原问题变形为
$$
\begin{array}{lrrrr}
\operatorname{minimize} & z \\
\text {subject to} & z & - & \mathbf{c}^T\mathbf{x} & = & 0\\
&&& \mathbf{Ax} & = & \mathbf{b}\\
&&& \mathbf{x} & \geq & \mathbf{0}
\end{array}
$$
选择基变量后，有
$$
\begin{array}{lrrrrrr}
\text {minimize} & z & & & & \\
\text {subject to} & z & - & \mathbf{c}_{B} \mathbf{x}_{B} & - & \mathbf{c}_{N} \mathbf{x}_{N} & = & 0 \\
& & & \mathbf{B x}_{B} & + & \mathbf{N} \mathbf{x}_{N} & = & \mathbf{b} \\
& & & \mathbf{x}_{B}& , & \mathbf{x}_{N} & \geq & \mathbf{0} .
\end{array}
$$
由于$\mathbf{B}$可逆，因此
$$
\mathbf{x}_{B}+\mathbf{B}^{-1} \mathbf{N} \mathbf{x}_{N}=\mathbf{B}^{-1} \mathbf{b}
$$
解得$\mathbf{x}_B$代入第一条约束中，可得
$$
z+\mathbf{0} \mathbf{x}_{B}+\left(\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{N}-\mathbf{c}_{N}\right) \mathbf{x}_{N}=\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{b}
$$

最终，得到了
$$
\begin{array}{lrrrrrr}
\text {minimize} & z & & & & \\
\text {subject to} & z & + & \mathbf{0} \mathbf{x}_{B} & + & (\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{N}-\mathbf{c}_{N}) \mathbf{x}_{N} & = & \mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{b} \\
& & & \mathbf{x}_{B} & + & \mathbf{B}^{-1} \mathbf{N} \mathbf{x}_{N} & = & \mathbf{B}^{-1} \mathbf{b} \\
& & & \mathbf{x}_{B}& , & \mathbf{x}_{N} & \geq & \mathbf{0} .
\end{array}
$$

将所有系数放入表格中，即有：

|                |      $z$     | $\mathbf{x}_B$ |                             $\mathbf{x}_N$                            |                         RHS                        |                  |
|:--------------:|:------------:|:--------------:|:---------------------------------------------------------------------:|:--------------------------------------------------:|:----------------:|
|       $z$      |       1      |  $\mathbf{0}$  | $\mathbf {c} _{B}  \mathbf {B}^{- 1 }  \mathbf {N}- \mathbf {c}_ {N}$ | $\mathbf {c}_{B}  \mathbf {B}^{- 1 }  \mathbf {b}$ |       Row 0      |
| $\mathbf{x}_B$ | $\mathbf{0}$ |  $\mathbf{I}$  |                   $\mathbf {B}^{- 1 }  \mathbf {N}$                   |          $\mathbf {B}^{- 1 }  \mathbf {b}$         | Rows 1 through m |

非基变量$\mathbf{x}_N = 0$时，可以在表格**RHS**一列得到此时的目标函数值和各基变量的值。每次迭代进基和出基时，只需要进行矩阵的初等行变换，将旧基下的该形式变换为新基下的该形式。这里，大部分教材都进行了详细举例，我们不妨略过。

最后，我们对进基和出基进行一些讨论。进出基需要满足下列三个要求：
- 进出基变换后，所有基变量还是一组基（即系数矩阵满秩）；
- 求得的基解可行（即$\mathbf{x}_B \geq 0$）；
- 目标函数值减少或不变化（对于最小化问题）。

若我们最初选择的基对应的基解是可行的，那么我们的进出基显然是满足这三个要求的。因为：
- $y_{rk}\neq 0$（注意到$\mathbf{B}\bar{\mathbf{a}_k} = \mathbf{a}_k$）；
- $r = \arg\min\limits_{i}\{\frac{\bar{b_i}}{y_{ik}}: y_{ik} > 0\}$；
- $\mathbf{c}_B^T\mathbf{\bar{a}}_{k} - {c}_k>0$。

那么，问题来了。第一，如何选择满足要求的初始基？第二，是否可以一次性选择多个进基变量和出基变量。对于第一个问题，我们将在下一章进行回答。对于第二个问题，难点在于第二个要求（即求得的基解仍然可行）往往难以验证。除了一些特殊问题外，一次性选择多个进基变量和出基变量往往是被禁止的。