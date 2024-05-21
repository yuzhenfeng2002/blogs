# 矩阵的标准形

## 相似矩阵

考虑关于$X$的方程$AX = \lambda X$，即：
$$
AX - \lambda X = AX - \lambda EX = (A - \lambda E)X = 0
$$
是一个齐次线性方程组，我们需要考察该齐次线性方程组是否有非零解。显然，当
$$
\mathrm{R}(A-\lambda E) < n \Leftrightarrow \det(A-\lambda E) = 0
$$
时，有非零解。$\det(A-\lambda E)$是一个$\lambda$的$n$阶多项式，$n$是方阵$A$的阶数。$\det(A-\lambda E) = 0$这个$n$阶多项式的零点在复数范围内恰有$n$个（代数学基本定理），这$n$个零点称为$A$的特征值，相应的非零解称为关于该特征值$\lambda$的特征向量。

**定理** 不同特征值对应的特征向量线性无关。

假设$A$恰好有$n$个不同的特征值$\lambda_i$，每个特征值对应一个特征向量$p_i$，根据特征值特征向量的定义有：
$$
Ap_i = \lambda_ip_i
$$
记$P = (p_1, p_2, \cdots, p_n)$，则有
$$
AP = P\cdot\mathrm{diag}(\lambda_1, \lambda_2, \cdots, \lambda_n) = P\Lambda
$$
由于特征向量$p_i$线性无关，则有
$$
A = P\Lambda P^{-1}
$$
$A$与对角阵$\Lambda$相似。

**例** $A^n$
$$
A^n = (P\Lambda P^{-1})^n = P\Lambda^n P^{-1}
$$

如果$A$有相同的特征值呢？例如矩阵$\begin{pmatrix}2 & 1 \\ 0 & 2 \end{pmatrix}$有两特征值$\lambda_1 = \lambda_2 = 2$，对应的无关的特征向量只有$p = c(1, 0)^T$。此时就无法找到前面$A = P\Lambda P^{-1}$中与$A$相似的$\Lambda$。因此，这一章的基本问题是对于给定矩阵$A$，寻找“相对好的”$J$使得$A$与$J$相似。

## 一元多项式

形如$\sum\limits_{k = 0}^{n} a_{k}\lambda^{k}, a_k \in \mathbb{K}$的表达式称为关于$\lambda$的一元多项式，记为$f(\lambda)$。若$a_n \neq 0$，则称$n$为该多项式的次数，记为$\deg(f(x))$；若$a_k = 0$，则称该多项式为零多项式，零多项式的次数可以定义为$-\infty$。数域$\mathbb{K}$上一切一元多项式构成的集合，记为$\mathbb{F}[\lambda]$。

多项式间的基本运算包括加法和乘法。其运算有如下性质：

**性质** $\deg(f(x) + g(x)) \leq \max(\deg(f(x)), \deg(g(x)))$

**性质** $\deg(f(x) g(x)) = \deg(f(x)) + \deg(g(x))$

其中加法具有逆运算（即减法），乘法没有逆运算（即除法），但是有带余除法运算。例如现有$f(x) = x^4 + x^3 + 1$，$g(x) = x - 2$，则有：
$$
f(x) = g(x)q(x) + r(x), q(x) = x^3 + 3x^2 + 6x + 12, r(x) = 25
$$
其中$q(x)$称为商式，$r(x)$称为余式。若$r(x) = 0$，则称$f(x)$可以被$g(x)$整除，此时$f(x)$称为$g(x)$的倍式，$g(x)$称为$f(x)$的因式。

**最大公因式** 若$d(x)$为$f(x)$和$g(x)$的公因式，且对一切$f(x)$和$g(x)$的公因式$\widetilde{d}(x)$，都有$d(x)$能被$\widetilde{d}(x)$整除，则称$d(x)$为$f(x)$和$g(x)$的最大公因式。这和数的最大公因式的定义基本是一致的。

最小公倍式可仿照最大公因式进行定义。

寻找一个数的最大公因数时，我们可以使用素因数分解，因为一定范围内的素数的个数是有限的。但是，我们求多项式的最大公因式时，素因数这个方法显然失效了——一定范围内的不可约多项式的个数是无限的。这个时候应该使用**辗转相除法**。

**证**
$$
    \begin{cases}
        f(x) = g(x)q(x) + r(x)\\
        g(x) = r(x)q_1(x) + r_1(x)\\
        r(x) = r_1(x)q_2(x) + r_2(x)\\
        \cdots\\
        r_{n-2}(x) = r_{n-1}(x)q_{n}(x) + r_{n}(x)\\
        r_{n-1}(x) = r_{n}(x)q_{n+1}(x)
    \end{cases}
$$
从下往上依次可证$r_n$可以被$r_{n-1}(x), r_{n-2}(x),\cdots, g(x), f(x)$整除；假设有$f(x)$和$g(x)$的公因式$\widetilde{d}(x)$，则从上往下依次可证$\widetilde{d}(x)$可以被$r(x), r_{1}(x),\cdots, r_{n}(x)$整除。

两多项式的众多最大公因式只相差一个常数因子（根据定义，$d(x)$和$\widetilde{d}(x)$能互相被整除），其中最高次数项系数为1的最大公因式记为$\gcd(f(x), g(x))$。最小公倍式即为$\frac{f(x)g(x)}{d(x)}$，记作$\mathrm{lcm}(f(x), g(x))$。

通过从上往下消去$r(x), r_1(x), \cdots, r_{n-1}(x)$可以如下性质：

**性质** $d(x) = u(x)f(x) + v(x)g(x), d(x) = \gcd(f(x), g(x))$

**例** $4x + 12y = 10$有无整数解？
$$
f(x) = 4, g(x) = 12, \gcd(f(x), g(x)) = 4, \mod{10, 4} = 2 \neq 0
$$

最大公因式与数域无关，但是因式分解跟数域有关。例如$x^2 + 1$在有理数域为不可约多项式，但在复数域便可以继续因式分解。在复数域因式分解，任何多项式都可以分解到一次因式。

**因式定理** 在复数域，$(x - a)$能够整除多项式$f(x)$。

## Lambda矩阵的标准形

**定义** 每个元素均为$\mathbb{F}[\lambda]$的元素的矩阵，称为$\lambda$矩阵。

$\lambda$矩阵的三类初等变换：

- 交换两行（列）的位置；
- 某行（列）乘一个非零**数**；
- 某行（列）加另一行（列）的**多项式**倍。

有限次的初等变换可以将一个数阵变换为
$
\begin{pmatrix}
    E_r & O\\
    O & O
\end{pmatrix}
$
变换的过程中，矩阵的秩没有发生改变。那么，初等变换对于$\lambda$矩阵呢？显然不行，例如对于
$
\begin{pmatrix}
    \lambda & 0\\
    0 & 0
\end{pmatrix}
$
，无法将其变成上述形式。但是在对$\lambda$矩阵进行初等变换时，$\lambda$矩阵所含元素的最大公因式没有发生改变。

**证** 对于前两种初等变换，显然；对第三种进行讨论：
$$
    d(x) = \gcd(f(x), g(x)) \Rightarrow d(x)|(f(x) + g(x)h(x))\\
    \text{Assume } \widetilde{d}(x)|(f(x) + g(x)h(x)), \widetilde{d}(x)|g(x), d(x)|\widetilde{d}(x)\\
    \text{Therefore } \widetilde{d}(x)|g(x), \widetilde{d}(x)|f(x)
$$
这与$d(x)$为$f(x)$和$g(x)$的最大公因式矛盾。

根据辗转相除的原理，$\lambda$矩阵$A_0$经过有限次的初等变换后，可以得到所有元素的最大公因式，不妨设出现在第1行第1列处，即
$
A_0 \rightarrow
\begin{pmatrix}
    d_1(\lambda) & O\\
    O & A_1
\end{pmatrix}
$。

$A_1$经过有限次的初等变换后，也可以得到剩余所有元素的最大公因式，不妨设出现在第2行第2列处，即
$
A_1 \rightarrow
\begin{pmatrix}
    d_1(\lambda) & 0 & O\\
    0 & d_2(\lambda) & O\\
    O & O & A_2
\end{pmatrix}
$。

于是可以由初等变换得到形如：
$$
\begin{pmatrix}
    \mathrm{diag}(d_1(\lambda), d_2(\lambda), \cdots, d_r(\lambda)) & O\\
    O & O
\end{pmatrix}
$$
的矩阵，其中$d_1(\lambda)|d_2(\lambda)|\cdots|d_r(\lambda)$。

类似初等变换保持同阶子行列式是否为零的特性不改变，我们有：

**定理** 对$\lambda$矩阵进行初等变换保持同阶子行列式的公因式不变，同阶子式的最大公因式称为该矩阵的行列式因子。

于是$D_1(\lambda) = d_1(\lambda)$为1阶行列式因子，$D_2(\lambda) = d_1(\lambda)d_2(\lambda)$为2阶行列式因子，以此类推。在实际问题中，有时行列式因子更容易得到，可以通过行列式因子间接求$d_i(\lambda)$。

**例**
已知$
A = \begin{pmatrix}
    1 & -1 & 0\\
    2 & 4 & -1\\
    0 & 0 & 3\\
\end{pmatrix}
$，求$\lambda E - A$的标准形。

**法一**
$$
\lambda E - A = \begin{pmatrix}
    \lambda - 1 & 1 & 0\\
    -2 & \lambda - 4 & 1\\
    0 & 0 & \lambda - 3
\end{pmatrix}
$$
观察到一阶行列式因子$d_1(\lambda) = 1$；又观察到子式
$\det \begin{pmatrix}
    1 & 0\\
    \lambda - 4 & 1
\end{pmatrix} = 1$，故二阶行列式因子$d_1(\lambda)d_2(\lambda) = 1$；三阶行列式因子$d_1(\lambda)d_2(\lambda)d_3(\lambda) = \det(\lambda E - A) = (\lambda - 3)^2(\lambda - 2)$。故$\lambda E - A$的标准形为$\mathrm{diag}(1, 1, (\lambda - 3)^2(\lambda - 2))$。

**法二**
$$
    \lambda E - A\rightarrow
    \begin{pmatrix}
        1 & \lambda - 1 & 0\\
        \lambda - 4 & -2 & 1\\
        0 & 0 & \lambda - 3
    \end{pmatrix}\rightarrow
    \begin{pmatrix}
        1 & \lambda - 1 & 0\\
        0 & -2 - (\lambda - 1)(\lambda - 4) & 1\\
        0 & 0 & \lambda - 3
    \end{pmatrix}\rightarrow\\
    \begin{pmatrix}
        1 & 0 & 0\\
        0 & 1 & -2 - (\lambda - 1)(\lambda - 4)\\
        0 & \lambda - 3 & 0
    \end{pmatrix}\rightarrow
    \begin{pmatrix}
        1 & 0 & 0\\
        0 & 1 & -2 - (\lambda - 1)(\lambda - 4)\\
        0 & 0 & (\lambda - 3)^2(\lambda - 2)
    \end{pmatrix}
$$

对于矩阵$\lambda E - A$而言，其标准形形如$\mathrm{diag}(d_1(\lambda), d_2(\lambda), \cdots, d_r(\lambda))$。$d_i(\lambda)$称为$A$的**不变因子**；将各不变因子在复数域上作因式分解，得到的一次因式连同其幂称为**初等因子**。对于分块矩阵，每一块的初等因子，恰是整个矩阵的初等因子。

## Jordan标准形

**定理**
$A$一定复相似于矩阵
$$
J = \mathrm{diag}(J_1, J_2, \cdots, J_s), 
J_k = \begin{pmatrix}
    \lambda_k & 1 & & \\
    & \ddots & \ddots & \\
    & & \lambda_k & 1\\
    & & & \lambda_k
\end{pmatrix}
$$

$J_k$称为Jordan块，$J = \mathrm{diag}(J_1, J_2, \cdots, J_s)$称为Jordan标准形。当$J_k$均为一阶时，即是对角阵。记$J_k$的阶数为$r_k$，矩阵$A$的初等因子分别为
$$
(\lambda - \lambda_1)^{r_1}, (\lambda - \lambda_2)^{r_2}, \cdots (\lambda - \lambda_s)^{r_s}
$$

**例**
已知$
A = \begin{pmatrix}
    1 & -1 & 0\\
    2 & 4 & -1\\
    0 & 0 & 3\\
\end{pmatrix}
$，求$A$的Jordan标准形。

$\lambda E - A$的标准形为$\mathrm{diag}(1, 1, (\lambda - 3)^2(\lambda - 2))$，初等因子为$(\lambda - 3)^2, (\lambda - 2)$，故$A$的Jordan标准形为
$$
\begin{pmatrix}
    3 & 1 & 0\\
    0 & 3 & 0\\
    0 & 0 & 2\\
\end{pmatrix}
$$

根据前面的经验，当矩阵$A$有$n$个不同的特征值时，
$$
A^n = (P\Lambda P^{-1})^n = P\Lambda^n P^{-1}
$$
引出Jordan标准形后，对于矩阵$A$，
$$
A^n = (PJ P^{-1})^n = PJ^n P^{-1}
$$

根据分块矩阵乘法，考察$J_k^n$：
$$
    J_k = \begin{pmatrix}
        \lambda_k & 1 & & \\
        & \ddots & \ddots & \\
        & & \lambda_k & 1\\
        & & & \lambda_k
    \end{pmatrix} = \lambda_kE + H\\
    J_k^n = (\lambda_kE + H)^n
$$
$H$矩阵乘法有特殊的性质，又由于$\lambda_kEH = H\lambda_kE$，所以$(\lambda_kE + H)^n$可以进行二项式展开。

那么，如何求**过渡矩阵**$P$呢？

**例**
已知$
A = \begin{pmatrix}
    1 & -1 & 0\\
    2 & 4 & -1\\
    0 & 0 & 3\\
\end{pmatrix}
$，求$A$的Jordan标准形及过渡矩阵$P$。

由上，$A$的Jordan标准形为

$$
J = \begin{pmatrix}
    3 & 1 & 0\\
    0 & 3 & 0\\
    0 & 0 & 2\\
\end{pmatrix}
$$

因此，
$$
A = PJP^{-1} \rightarrow AP = PA
$$
将矩阵$P$分块，则有$P = (p_1, p_2, p_3)$，于是，
$$
    AP = (Ap_1, Ap_2, Ap_3)\\
    PJ = (3p_1, p_1+3p_2, 2p_3)\\
    AP = PJ \rightarrow 
    \begin{cases}
        Ap_1 = 3p_1\\
        Ap_2 = p_1 + 3p_2\\
        Ap_3 = 2p_3
    \end{cases}
$$
可依次求得$P$的一个解为
$$
P = \begin{pmatrix}
    1& -1& 1\\
    -2& 1& -1\\
    0& 1& 0
\end{pmatrix}
$$

**例**
已知$
A = \begin{pmatrix}
    -1 & 2 & -2\\
    -1 & 2 & -1\\
    1 & -1 & 2
\end{pmatrix}
$，求$A$的Jordan标准形及过渡矩阵$P$。

$A$的Jordan标准形为
$$
J = \begin{pmatrix}
    1 & 0 & 0\\
    0 & 1 & 1\\
    0 & 0 & 1
\end{pmatrix}
$$

因此，
$$
    AP = (Ap_1, Ap_2, Ap_3)\\
    PJ = (p_1, p_2, p_2 + p_3)\\
    AP = PJ \rightarrow 
    \begin{cases}
        Ap_1 = p_1\\
        Ap_2 = p_2\\
        Ap_3 = p_2 + p_3
    \end{cases}
$$
由于$Ap_1 = p_1, Ap_2 = p_2$，所以$(A - E)p = 0$有两个线性无关解，如何选择$p_1$和$p_2$呢？$p_2$应使得方程$Ap_3 = p_2 + p_3$有解。

$$
    (A - E)p = 0 \rightarrow
    p = c_1(1, 1, 0)^T + c_2(1, 0, -1)^T\\
    (A - E)p_3 = (c_1 + c_2, c_1, -c_2)\\
    \begin{pmatrix}
        -2 & 2 & -2 & c_1 + c_2\\
        -1 & 1 & -1 & c_1\\
        1 & -1 & 1 & -c_2
    \end{pmatrix} \rightarrow
    \begin{pmatrix}
        0 & 0 & 0 & c_1 - c_2\\
        0 & 0 & 0 & c_1 - c_2\\
        1 & -1 & 1 & -c_2
    \end{pmatrix}
$$

因此，要$(A - E)p_3 = (c_1 + c_2, c_1, -c_2)$有解，需要$c_1 = c_2$。不妨取$c_1 = c_2 = 1$，则
$$
p_3 = (0, 1, 0)^T, p_2 = (2, 1, -1)^T
$$
$p_1$可任取与$p_2$线性无关的向量$(1, 1, 0)$。

因此求得$P$的一个解为
$$
P = \begin{pmatrix}
    1& 2& 0\\
    1& 1& 1\\
    0& -1& 0
\end{pmatrix}
$$

## 矩阵函数

**定义** $f(A)$
$$
f(A) = a_nA^n + a_{n-1}A^{n-1} + \cdots + a_0E
$$
则，
$$
f(A) = P(a_nJ^n + a_{n-1}J^{n-1} + \cdots + a_0E)P^{-1}
= Pf(J)P^{-1}
$$

下面研究$f(J)$的算法。

**法一**

$$
    J = \mathrm{diag}(J_1, J_2, \cdots, J_s)\\
    J^k = \mathrm{diag}(J_1^k, J_2^k, \cdots, J_s^k)\\
    J_0^k = (\lambda_0E + H)^k = \lambda_0^kE + k\lambda_0^{k-1}H + \frac{k(k-1)}{2}\lambda_0^{k-2}H^2 + \cdots\\
    = \begin{pmatrix}
        \lambda_0^k & k\lambda_0^{k-1} & \frac{k(k-1)}{2}\lambda_0^{k-2} & \cdots \\
        & \ddots & \ddots & \ddots\\
        % &  & \ddots & \ddots\\
        % &  &  & \ddots
    \end{pmatrix}
    = \begin{pmatrix}
        \lambda_0^k & (\lambda_0^k)^\prime & \frac{(\lambda_0^k)^{(2)}}{2!} & \cdots \\
        & \ddots & \ddots & \ddots\\
        % &  & \ddots & \ddots\\
        % &  &  & \ddots
    \end{pmatrix}\\
    f(J_0) = \begin{pmatrix}
        f(\lambda_0) & f(\lambda_0)^\prime & \frac{f(\lambda_0)}{2!} & \cdots \\
        & \ddots & \ddots & \ddots\\
        % &  & \ddots & \ddots\\
        % &  &  & \ddots
    \end{pmatrix}\\
    f(J) = \mathrm{diag}(f(J_1), f(J_2), \cdots, f(J_s))
$$

在这种方法下，过渡矩阵$P$可能较难求解。因此，介绍不需要求解过渡矩阵的法二。

**法二**

寻找**零化多项式**$g(x)$，使得$g(A) = 0$，求$f(x)$关于$g(x)$的余式$r(x)$，则$f(A) = r(A)$。

**Hamilton-Cayley定理**
$$ g(\lambda) = \det(\lambda E - A), g(A) = 0 $$

也可以取次数更低的$g(\lambda) = d_n(\lambda)$，但不存在次数比$g(\lambda) = d_n(\lambda)$更低的零化多项式，因此称其为**最小多项式**。

**解释**
对于每一个Jordan块，有$H = J_0 - \lambda_0 E$。由于$H^r = (J_0 - \lambda_0 E)^r = 0$，因此$(x - \lambda_0)^r$是矩阵$J_0$的零化多项式，$(x - \lambda_0)^r$是矩阵$J_0$的初等因子，$\det(\lambda E - A)$是所有初等因子的乘积。

此时，只需要计算比$f(A)$多项式次数小的多项式$r(A)$即可。

形式上定义矩阵的幂级数，即为
$$
f(A) = a_0E + a_1A + \cdots + a_nA^n + \cdots = Pf(J)P^{-1}
$$

求矩阵的幂级数时，可以使用上述求解矩阵多项式的方法一。由此可以得到矩阵的幂级数收敛的充分必要条件，即矩阵的特征值都落在幂级数的收敛域中。

**例**
$$
e^A = E + A + \frac{A^2}{2!} + \cdots
$$

**例**
已知
$$
    A = \begin{pmatrix}
        -1 & 2 & -2\\
        -1 & 2 & -1\\
        1 & -1 & 2
    \end{pmatrix},
    J = \begin{pmatrix}
        1 & 0 & 0\\
        0 & 1 & 1\\
        0 & 0 & 1
    \end{pmatrix},
    P = \begin{pmatrix}
        1& 2& -1\\
        1& 1& 0\\
        0& -1& 0
    \end{pmatrix}
$$
求$e^A$。

$$
    e^A = Pe^JP^{-1} = P\begin{pmatrix}
        e & 0 & 0\\
        0 & e & e\\
        0 & 0 & e
    \end{pmatrix} = P(eJ)P^{-1} = eA
$$

仿照上述求解矩阵多项式的方法二，是否可以利用其来求解矩阵的幂级数呢？

记$g(x) = (x-x_1)^{\alpha_1}(x-x_2)^{\alpha_2}\cdots (x-x_s)^{\alpha_s}$，可将$f(x)$在$x_k$处展开，于是
$$
f(x) \equiv r_k(x)(\mod (x - x_k)^{\alpha_k})
$$

**中国剩余定理**
存在唯一的$r(x)$，次数小于$\deg g(x)$，使得
$$
r(J_k) = r_k(J_k)(\mod (x - x_k)^{\alpha_k})
$$
此时，$r(J_k) = f(J_k)$，因此$f(A) = r(A)$。

**例**
已知$
A = \begin{pmatrix}
    1 & -1 & 0\\
    2 & 4 & -1\\
    0 & 0 & 3\\
\end{pmatrix}
$，求$e^A$。

$$
    \lambda E - A\rightarrow
    \begin{pmatrix}
        1 & 0 & 0\\
        0 & 1 & -2 - (\lambda - 1)(\lambda - 4)\\
        0 & 0 & (\lambda - 3)^2(\lambda - 2)
    \end{pmatrix}\\
    e^{\lambda} \equiv e^3 + e^3(\lambda - 3)(\mod (\lambda - 3)^2)\\
    e^{\lambda} \equiv e^2 (\mod (\lambda - 2))
$$
利用待定系数法，设$r(\lambda) = e^3 + e^3(\lambda - 3) + a(\lambda - 3)^2$，则由于$r(\lambda) = (\lambda - 2)q(\lambda) + \widetilde{r}(\lambda)$
$$
    r(\lambda) \equiv \lambda(2)(\mod (\lambda - 2))\\
    r(2) = e^3 - e^3 + a = e^2\rightarrow a = e^2
$$

**定理**
$A$可对角化的充分必要条件是$A$的最小多项式$d_n(\lambda)$无重根。

**定理**
$f(\lambda)$无重根的充分必要条件是$\gcd (f(\lambda), f^\prime(\lambda)) = 1$。

**解释**
$$
    f(x) = (x - x_1)^n g(x)\\
    f^\prime(x) = (x - x_1)^{n-1}(ng(x) + (x - x_1)g^\prime(x))
$$

## 矩阵函数与微分方程组的解

**常系数线性常微分方程组**
设$x_1(x), x_2(x), \cdots, x_n(x)$未知，
$$
\begin{cases}
    x_1'(t) = a_{11}x_1(t) + a_{12}x_2(t) + \cdots + x_{1n}x_n(t)\\
    x_2'(t) = a_{21}x_1(t) + a_{22}x_2(t) + \cdots + x_{2n}x_n(t)\\
    \cdots\\
    x_n'(t) = a_{n1}x_1(t) + a_{n2}x_2(t) + \cdots + x_{nn}x_n(t)
\end{cases}
$$

令$X(t) = (x_1(t), x_2(t), \cdots, x_n(t))^T, A=(a_{ij})$，则上述方程组即为
$$
X'(t) = AX(t)
$$

**验证一** $X(t) = e^{At}C$为解
$$
    e^{At} = \sum\limits_{n=0}^{\infty} \frac{A^n t^n}{n!}\\
    \frac{\mathrm{d}e^{At}}{\mathrm{d}t}
    = \sum\limits_{n=0}^{\infty} \frac{A^n nt^{n-1}}{n!}
    = \sum\limits_{n=0}^{\infty} \frac{A^{n+1} t^n}{n!}\\
    = Ae^{At} = e^{At}A\\
    \frac{\mathrm{d}(e^{At}C)}{\mathrm{d}t} = Ae^{At}C
$$

**验证二**
$$
    X'(t) = AX(t)\rightarrow e^{-At}X'(t) = e^{-At}AX(t)\\
    \rightarrow e^{-At}X'(t) - e^{-At}AX(t) = \frac{\mathrm{d}e^{-At}X(t)}{\mathrm{d}t} = 0\\
    e^{-At}X(t) \equiv C, C = X(0)
$$

**例**
已知$
A = \begin{pmatrix}
    3 & 0 & 8\\
    3 & -1 & 6\\
    -2 & 0 & -5\\
\end{pmatrix}
$，求$X'(t) = AX(t)$的通解。

通解为$X(t) = e^{At}C$，以下求$e^{At}$。设$f(\lambda) = e^{\lambda t}$，则$f(A) = e^{At}$。
$$
    \lambda E - A \rightarrow \mathrm{diag}(1, \lambda + 1, (\lambda + 1)^2)\\
    e^{\lambda t} \equiv e^{-t} + te^{-t}(\lambda + 1)(\mod{(\lambda + 1)^2})\\
    e^{At} = te^{-t}A + (e^{-t} + e^{-t}t)E
$$

## 习题选讲

**例1**
$A = \begin{pmatrix}
    0 & 1 &  & & \\
    0 & 0 & 1 & & \\
     &  &  & \ddots \\
     &  &  &  & 1\\
    -a_0 & -a_1 & \cdots & -a_{n-2} & -a_{n-1}\\
\end{pmatrix}, a_i\in \mathbb{R}$，求$\lambda E - A$的标准形。
$$
    \lambda E - A = \begin{pmatrix}
        \lambda & -1 &  & & \\
        0 & \lambda & -1 & & \\
         &  & \ddots & \ddots \\
         &  &  & \lambda & -1\\
        a_0 & a_1 & \cdots & a_{n-2} & \lambda+a_{n-1}\\
    \end{pmatrix}
$$
$\lambda E - A$存在一个$n - 1$阶子式为非零常数，故
$$
    \lambda E - A \rightarrow \mathrm{diag}(1, 1, \cdots, 1, \det(\lambda E - A))\\
    \det(\lambda E - A) = \lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_0
$$
进一步考察$A$矩阵，有
$$
    X'(t) = AX(t)\\
    \Leftrightarrow
    \begin{cases}
        x_1'(t) = x_2(t)\\
        x_2'(t) = x_3(t)\\
        \cdots\\
        x_{n-1}'(t) = x_n(t)\\
        x_n'(t) = -a_0x_1(t) - a_1x_2(t) - \cdots - a_{n-1}x_n(t)
    \end{cases}\\
    \Leftrightarrow
    y^{(n)}(t) + a_{n-1}y^{(n-1)}(t) + \cdots + a_0y(t) = 0
$$

**例2** $A^2 = E$，证明$A$相似于对角阵。

**证** 由于$A^2 - E = 0$，有$\lambda ^2 - 1$为$A$的零化多项式，故$m_A(\lambda)|(\lambda - 1)(\lambda + 1)$，故$m_A(\lambda)$无重根，证毕。

**例3** $A^2 = E$，求$A$所有可能的Jordan标准形。

由于$m_A(\lambda)|(\lambda - 1)(\lambda + 1)$，故$J$可能的情况有：
$$
    \begin{pmatrix}
        1 &  & \\
         & \lambda\pm1  & \\
         &  & (\lambda - 1)(\lambda + 1)\\
    \end{pmatrix}, 
    \begin{pmatrix}
        \lambda \pm 1 &  & \\
         & \lambda \pm 1  & \\
         &  & \lambda \pm 1\\
    \end{pmatrix}
$$