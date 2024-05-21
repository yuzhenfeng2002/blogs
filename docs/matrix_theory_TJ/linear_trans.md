# 线性空间与线性变换

## 线性空间

记$\mathbb{R}^2$为一切二维向量全体，对加法和数乘两种运算封闭；同样矩阵$\mathbb{F}^{m\times n}$也对和数乘两种运算封闭。它们满足8条基本性质：

- $u + v = v + u$
- $(v_1 + v_2) + v_3 = v_1 + (v_2 + v_3)$
- $0 + v = v$
- $v + (-v) = 0$
- $1v = v$
- $k(lv) = (kl)v$
- $(k + l)v = kv + lv$
- $k(u + v) = ku + kv$


**定义**
设$\mathbb{F}$是一个数域，$V$是一个集合，$V$上有两种运算，称为加法和数乘，满足以上8条性质，称$V$为数域$\mathbb{F}$上的线性空间，$V$中的元素称为向量。

**例** $\mathbb{F}[x]$

**例** $[a, b]$上的连续函数全体的集合，记作$C([a, b])$

**例** $V = \mathbb{R}^+, \mathbb{F} = \mathbb{R}$。定义$a\oplus b = ab, k\otimes a = a^k$

**例** 证明有且仅有一个$-v$，满足$v + (-v) = 0$

$$
    v + u_1 = 0, v + u_2 = 0\\
    \Rightarrow u_1 = u_1 + 0 = u_1 + (v + u_2) = (u_1 + v) + u_2 = u_2
$$

**例** 证明$-v = (-1)v$

$$
    (-1)v + v = (-1)v + 1v = (-1 + 1)v = 0v\\
    0v + v = 0v + 1v + (0+1)v = 1v = v\\
    \Rightarrow(0v + v) + (-v) = v + (-v) = 0\\
    \Rightarrow LHS = 0v + (v + (-v)) = 0v + 0 = 0v = RHS = 0
$$

定义一个线性空间需要规定$(V, \mathbb{F}, \oplus, \otimes)$。

## 线性空间的维数、基与坐标

**线性相关** 设$v_1, v_2, \cdots, v_n \in V$，若存在$k_1, k_2, \cdots, k_n \in \mathbb{F}$不全为0，且
$$
\sum\limits_{i = 1}^n k_iv_i = 0
$$
则称$v_1, v_2, \cdots, v_n$线性相关，否则称线性无关。

**有限维线性空间** 存在有限个极大线性无关组的线性空间。

设$V$为有限维线性空间，$v_1, v_2, \cdots, v_n$为极大线性无关组。任取$v \in V$，于是存在$k_1, k_2, \cdots, k_n, k \neq 0\in \mathbb{F}$，使得
$$
    kv = \sum\limits_{i = 1}^n k_iv_i\\
    \Rightarrow v = \sum\limits_{i=1}^n \frac{k_i}{k}v_i
$$

在此基础上，可以得到，任取$v \in V$，存在唯一的$c_1, c_2, \cdots, c_n \in \mathbb{F}$，使得
$$
    v = \sum\limits_{i = 1}^n c_iv_i
$$
$v_1, v_2, \cdots, v_n$称为$V$的一组**基**，$c_1, c_2, \cdots, c_n$为$v$在基$v_1, v_2, \cdots, v_n$下的**坐标**，记为$\mathrm{crd}(v: v_1, v_2, \cdots, v_n)$。

**例** $\mathbb{R}^n$

若以$e_1, e_2, \cdots, e_n$为基，则向量$(x_1, x_2, \cdots, x_n)$在该组基下的坐标仍然为$(x_1, x_2, \cdots, x_n)$。

**例** $\mathbb{R}_3[x]$

若以$1, x, x^2, x^3$为基，则向量$f(x) = 2x + x^3$在该组基下的坐标为$(0, 2, 0, 1)$。

**例** $\mathbb{R}^{2\times 2}$
若以
$$
\begin{pmatrix}1 & 0\\0 & 0\end{pmatrix},
\begin{pmatrix}0 & 1\\0 & 0\end{pmatrix},
\begin{pmatrix}0 & 0\\1 & 0\end{pmatrix},
\begin{pmatrix}0 & 0\\0 & 1\end{pmatrix},$$
为基，则向量$\begin{pmatrix}1 & 2\\3 & 4\end{pmatrix}$在该组基下的坐标为$(1, 2, 3, 4)$。

**例** $\mathbb{C}$为数域$\mathbb{C}$上的一维线性空间。

**例** $\mathbb{C}$为数域$\mathbb{R}$上的二维线性空间。

下面考察在两组不同的基下同一向量的坐标。即
$$
    (x_1, x_2, \cdots, x_n) = \mathrm{crd}(u: v_1, v_2, \cdots, v_n)\\
    (y_1, y_2, \cdots, y_m) = \mathrm{crd}(u: w_1, w_2, \cdots, w_m)
$$

**证明** $m = n$

用基$v_1, v_2, \cdots, v_n$分别表示$y_1, y_2, \cdots, y_m$，由于$y_1, y_2, \cdots, y_m$线性无关，所以$m \leq n$；反之，也有$n \leq m$，因此$m = n$。

设$w_j = \sum\limits_{i = 1}^n t_{ij}v_i$，则
$$
u = \sum\limits_{j = 1}^n\sum\limits_{i = 1}^n y_jt_{ij}v_i
= \sum\limits_{i = 1}^n(\sum\limits_{j = 1}^n y_jt_{ij})v_i
$$
因此，$x_i = \sum\limits_{j = 1}^n y_jt_{ij}$。

记$T = (t_{ij})_{n\times n} = (\mathrm{crd}(w_1: v), \cdots, \mathrm{crd}(w_n: v))$，于是有
$$
\mathrm{crd}(u: v_1, v_2, \cdots, v_n) = T\mathrm{crd}(u: w_1, w_2, \cdots, w_m)
$$
$T$称为$v$到$w$的**过渡矩阵**。

**例**
$$
    \alpha_1 = (1, 2, 1)^T, \alpha_2 = (2, 3, 3)^T, \alpha_3 = (3, 7, 10)^T\\
    \beta_1 = (3, 1, 4)^T, \beta_2 = (5, 2, 1)^T, \beta_3 = (1, 1, -6)^T
$$
首先，验证这两组向量组都是$\mathbb{R}^3$的一组基，即验证下面两矩阵的秩为$3$，利用初等变换将其变换为行阶梯矩阵即可。
$$
    \begin{pmatrix}
        1 & 2 & 3\\
        2 & 3 & 7\\
        1 & 3 & 10
    \end{pmatrix},
    \begin{pmatrix}
        3 & 5 & 1\\
        1 & 2 & 1\\
        4 & 1 & -6
    \end{pmatrix}
$$

再求从$\alpha$到$\beta$的过渡矩阵，即求解下列线性方程组。
$$
    \begin{pmatrix}
        1 & 2 & 3\\
        2 & 3 & 7\\
        1 & 3 & 10
    \end{pmatrix}X =
    \begin{pmatrix}
        3 & 5 & 1\\
        1 & 2 & 1\\
        4 & 1 & -6
    \end{pmatrix}
$$

记从$\alpha$到$\beta$的过渡矩阵为$T_{\alpha\beta}$，则有：

- $T_{\alpha\alpha} = E$；
- $T_{\alpha\beta} = A^{-1}B, T_{\beta\alpha} = B^{-1}A$，因此$T_{\alpha\beta}T_{\beta\alpha} = E$；
- $T_{\alpha\beta}T_{\beta\gamma} = T_{\alpha\gamma}$。


再次考察上述例题，$T_{\alpha\beta} = T_{\alpha e}T_{e\beta}$，其中$e$为标准基。$T_{e\beta}$可以轻易得到，又有$T_{\alpha e}T_{e\alpha} = E$，于是$T_{\alpha\beta} = A^{-1}B$。

## 子空间的运算
**定义** 假设$V$为$\mathbb{F}$上的线性空间，$W \subset V, W \neq \varnothing$，且$W$**关于$V$的运算构成**线性空间，则称$W$为$V$的子空间。

任意线性空间$V$均有两个子空间$V$与${0}$。

**例** $\mathbb{C}, \mathbb{R}$不是$\mathbb{C}, \mathbb{F}$的子空间（数域不同，数乘定义也就不同）。

**例** 过原点的直线或平面是三维空间的子空间。

**例** 若$v_1, v_2, \cdots, v_n \in V$，由$v_1, v_2, \cdots, v_n$生成的子空间$\mathrm{span}\{\sum\limits_{i = 0}^r k_iv_i, k_i \in \mathbb{F} \}$。

设$V_1, V_2$是$V$的两个子空间，则：

- $V_1\cap V_2$为子空间；
- $V_1 + V_2 = \{w_1 + w_2| w_1\in V_1, w_2\in V_2\}$；


**定理** $\dim (V_1 + V_2) + \dim (V_1\cap V_2) = \dim V_1 + \dim V_2$

**证明** 设$V_1\cap V_2$的一组基为$w_1, w_2, \cdots, w_r$；对于子空间$V_1$，进行扩充$w_{r+1}, \cdots, w_m$；同样对于子空间$V_2$，在$w_1, w_2, \cdots, w_r$基础上进行扩充$w_{m+1}, \cdots, w_n$。于是只要说明$V_1 + V_2$中的元素可以由$w_1, w_2, \cdots, w_n$线性表示，同时$w_1, w_2, \cdots, w_n$线性无关即可。

$V_1 + V_2$中的元素可以由$w_1, w_2, \cdots, w_n$线性表示是显然的。假设$w_1, w_2, \cdots, w_n$线性相关，则有
$$
\sum\limits_{i = 1}^r x_iw_i + \sum\limits_{i = 1}^{m - r} y_iw_{r + i} + \sum\limits_{i = 1}^{n - m} z_iw_{m + i} = 0
$$
也就是
$$
\sum\limits_{i = 1}^r x_iw_i + \sum\limits_{i = 1}^{m - r} y_iw_{r + i} = - \sum\limits_{i = 1}^{n - m} z_iw_{m + i}
$$
$\sum\limits_{i = 1}^r x_iw_i + \sum\limits_{i = 1}^{m - r} y_iw_{r + i} \in V_1, - \sum\limits_{i = 1}^{n - m} z_iw_{m + i} \in V_2$，所以$\sum\limits_{i = 1}^r x_iw_i + \sum\limits_{i = 1}^{m - r} y_iw_{r + i} \in V_1\cap V_2, y_i = 0$。同理，$z_i = 0$。因此$\sum\limits_{i = 1}^r x_iw_i = 0, x_i = 0$。得证。

**定义** 若子空间$V_1 \cap V_2 = \{0\}$，则$V_1 + V_2$称为直和，记作$V_1 \oplus V_2$。

**定理** 以下结论相互等价：

- $V_1 + V_2$为直和；
- 对一切$v \in V_1 + V_2$，均存在唯一的$v_1 \in V_1, v_2 \in V_2$使得$v = v_1 + v_2$；
- 若$v_1 + v_2 = 0, v_1 \in V_1, v_2 \in V_2$，则$v_1 = v_2 = 0$；
- 若$v_1, \cdots, v_k$为$V_1$的基，$v_{k+1}, \cdots, v_n$为$V_2$的基，则$v_1, \cdots, v_n$为$V_1 + V_2$的基；
- $\dim V_1 + \dim V_2 = \dim (V_1 + V_2)$。


**证明** 依次证明下列命题：

(1)推出(2)：设$v = v_1 + v_2 = v_1' + v_2', v_1, v_1'\in V_1, v_2, v_2'\in V_2$，则$v_1 - v_1' = v_2 - v_2'$。又由于$v_1 - v_1'\in V_1, v_2 - v_2'\in V_2$，所以$v_1 - v_1' = v_2 - v_2' \in V_1\cap V_2$。由于$V_1 + V_2$为直和，$V_1\cap V_2 = \{0\}$。

(2)推出(3)：(3)是(2)的特例。

(3)推出(4)：先证$v_1, \cdots, v_n$线性无关。设$\sum\limits_{i=1}^n l_iv_i = \sum\limits_{i=1}^k l_iv_i + \sum\limits_{i=k+1}^n l_iv_i = 0$，因此$\sum\limits_{i=1}^k l_iv_i = \sum\limits_{i=k+1}^n l_iv_i = 0$，故$l_i = 0$，又由于$\forall v \in V_1 + V_2$，有$v = v_1 + v_2, v_1 \in V_1, v_2 \in V_2$，故得证。

由(4)，(5)显然。

(5)推出(1)：根据定理$\dim (V_1 + V_2) + \dim (V_1\cap V_2) = \dim V_1 + \dim V_2$，有$\dim (V_1\cap V_2) = 0$，因此$V_1\cap V_2 = \emptyset$。

可将该定理扩展到$n$个子空间的直和。

## 线性变换

首先，考察正比例函数$x\mapsto y, f: \mathbb{R}\to \mathbb{R}$，满足
$$
f(x + y) = f(x) + f(y), f(lx) = lf(x)
$$

类似地，考虑线性空间到其自身的映射。

**定义** 若$\mathcal{A}: V \mapsto V$满足（1）$\forall v_1, v_2 \in V, \mathcal{A}(v_1 + v_2) = \mathcal{A}v_1 + \mathcal{A}v_2$，（2）$\forall v \in V, \mathcal{A}(kv) = k\mathcal{A}v$，则$\mathcal{A}$称为$V$上的线性变换。$V$上一切线性变换的全体，记作$\mathrm{End}(V)$。

以下例子均为线性变换。

**例** $V = \mathbb{R}^2$，（1）$v\mapsto v$逆时针旋转$\theta$角，（2）$v\mapsto v$关于x轴的镜像。

**例** $V = \mathbb{R}^{n\times 1}, A \in \mathbb{R}^{n\times n}, X\mapsto AX$。

**例** $V = \mathbb{F}[x], f(x)\mapsto f^\prime(x)$。

**例** $V = C([a, b]), f(x)\mapsto \int_{a}^{b} f(x) \mathrm{d}x$。

**例** 恒同变换：$v\mapsto v$；零变换：$v\mapsto 0$。

记$\mathcal{A}: V\to V, v\mapsto \mathcal{A}v, \mathrm{crd}v\mapsto \mathrm{crd}\mathcal{A}v$，先考虑基$v_1, \cdots, v_n$在该映射下的坐标变化。任取一向量$v$，假设$\mathrm{crd}v = (x_1, \cdots, x_n)^T$，则有
$$
    \mathrm{crd}\mathcal{A}v = \mathrm{crd}\mathcal{A}(x_1v_1 + \cdots + x_nv_n) = \mathrm{crd}(x_1\mathcal{A}v_1 + \cdots + x_n\mathcal{A}v_n) \\
    = x_1\mathrm{crd}\mathcal{A}v_1 + \cdots + x_n\mathrm{crd}\mathcal{A}v_n
    = (\mathrm{crd}\mathcal{A}v_1, \cdots, \mathrm{crd}\mathcal{A}v_n)(x_1, \cdots, x_n)^T
$$
记$A = (\mathrm{crd}\mathcal{A}v_1, \cdots, \mathrm{crd}\mathcal{A}v_n), \mathrm{crd}\mathcal{A}x = A\mathrm{crd}x$。$A$称为线性变换$\mathcal{A}$在基$v_1, \cdots, v_n$下的矩阵。

下面考察上述例子中的线性变换的矩阵。

**例** 恒同变换的矩阵为单位阵；零变换的矩阵为零矩阵。

**例** $V = \mathbb{R}^{n\times 1}, A \in \mathbb{R}^{n\times n}, X\mapsto AX$。该线性变换在基$e_1, \cdots, e_n$下的矩阵为$A$。

**例** $V = \mathbb{R}^2$，$v\mapsto v$逆时针旋转$\theta$角。该线性变换在基$e_1, e_2$下的矩阵为$\begin{pmatrix}
    \cos \theta & - \sin \theta\\
    \sin \theta & \cos \theta
\end{pmatrix}$。

**例** $V = \mathbb{R}^2$，$v\mapsto v$关于x轴的镜像。该线性变换在基$e_1, e_2$下的矩阵为$\begin{pmatrix}
    1 & 0\\
    0 & -1
\end{pmatrix}$。

**例** $\mathbb{F}_3[x], f(x)\mapsto f^\prime(x)$。该线性变换在基$1, x, x^2$下的矩阵为$\begin{pmatrix}
    0 & 1 & 0\\
    0 & 0 & 2\\
    0 & 0 & 0
\end{pmatrix}$。

下面考虑同一线性空间的两组基$v_1, \cdots, v_n$与$w_1, \cdots, w_n$，因此有
$$
\mathrm{crd}(u: v_1, v_2, \cdots, v_n) = T\mathrm{crd}(u: w_1, w_2, \cdots, w_n)
$$
$T$为过渡矩阵。$u = \mathcal{A}x$时，有$A_v\mathrm{crd}(x, v) = TA_w\mathrm{crd}(x, w)$。因此
$$
    A_v T\mathrm{crd}(x, w) = TA_w\mathrm{crd}(x, w)\\
    \Rightarrow T^{-1}A_v T\mathrm{crd}(x, w) = A_w\mathrm{crd}(x, w), \forall x
    \Rightarrow A_v = TA_wT^{-1}
$$

可以得到$A_v$相似于$A_w$。这就将该章与上一章相联系。

**线性变换的特征值**
对$\lambda \in \mathbb{F}$，存在$v \in V, v\neq 0$，使得$\mathcal{A}v = \lambda v$，则称$\lambda$为$\mathcal{A}$的特征值，$v$是关于$\lambda$的特征向量。

## 线性变换与子空间

由线性变换$\mathcal{A}$可以得到两个子空间，分别为：
- 核空间：$\ker \mathcal{A} = \{v| \mathcal{A}v = 0\}$
- 像空间：$\Im \mathcal{A} = \{w| w = \mathcal{A}v\}$

在一组基下，假设该线性变换的矩阵为$A$，核空间可以通过解线性方程组$Ax = 0$解得，像空间即矩阵$A$的列向量生成的线性空间。

**定义** 记$\mathcal{A}: V\to V, W\subseteq V$，若$\forall w\in W, \mathcal{A}w\in W$，则称$W$为$V$关于$\mathcal{A}$的不变子空间。

根据定义，可以看到$\ker \mathcal{A}$是$\mathcal{A}$的不变子空间。零空间与线性空间本身是一个不变子空间。

**定理** 记$V = V_1 \oplus V_2$，$\mathcal{A}$为$V$上的线性变换，且$V_1, V_2$为$\mathcal{A}$的不变子空间。若$w_1, \cdots, w_k$为$V_1$的基，$w_{k+1}, \cdots, w_n$为$V_2$的基，则$\mathcal{A}$在$w_1, \cdots, w_n$下的矩阵为分块对角阵。

**证明** 任取$w_j\in V_1$，则有$\mathcal{A}w_j\in V_1$，所以$\mathcal{A}w_j = \sum\limits_{i = 1}^k l_iw_i$，即$\mathcal{A}w_j$在基$w_1, \cdots, w_n$下的坐标为$(l_1, \cdots, l_k, 0, \cdots, 0)$。

## 习题选讲

考察集合$\mathbb{S}^{2}$。

**例1** 证明$\mathbb{S}^{2}$关于矩阵的加法和数乘构成线性空间。

由于$\mathbb{S}^{2}\subset \mathbb{R}^{2\times 2}$，而$\mathbb{R}^{2\times 2}$构成线性空间，因此只需验证$\mathbb{S}^{2}$关于加法和数乘封闭，这显然。

**例2** 求$\mathbb{S}^{2}$的一组基。

$\begin{pmatrix} 1 & 0\\ 0 & 0\end{pmatrix}, \begin{pmatrix} 0 & 0\\ 0 & 1\end{pmatrix}, \begin{pmatrix} 0 & 1\\ 1 & 0\end{pmatrix}$即为$\mathbb{S}^{2}$的一组基，接着证明该命题即可。

**例3** 在上一题的基下，求$\begin{pmatrix} 1 & 2\\ 2 & 0\end{pmatrix}$的坐标。

通过观察，坐标为$(1, 0, 2)^T$。

**例4** 考虑$\mathbb{R}^{2\times 2}$上的变换$\mathcal{T}: A\mapsto A^T$，证明$T$是$\mathbb{R}^{2\times 2}$上的线性变换。

$\mathcal{T}(A + B) = (A + B)^T = A^T + B^T = \mathcal{T}A + \mathcal{T}B, \mathcal{T}(kA) = (kA)^T = kA^T = k\mathcal{T}A$

**例5** 选$E_{11}, E_{12}, E_{21}, E_{22}$为$\mathbb{R}^{2\times 2}$的基，求$T$对应的矩阵。

由于$\mathcal{T}E_{11} = E_{11}, \mathcal{T}E_{12} = E_{21}, \mathcal{T}E_{21} = E_{12}, \mathcal{T}E_{22} = E_{22}$，于是
$$
T = \begin{pmatrix}
    1 & 0 & 0 & 0\\
    0 & 0 & 1 & 0\\
    0 & 1 & 0 & 0\\
    0 & 0 & 0 & 1\\
\end{pmatrix}
$$

**例6** 求$\mathcal{T}$的特征值。

即求$T$的特征值，即解特征多项式$|\lambda E - T| = (\lambda - 1)^3(\lambda + 1) = 0$。解得$\lambda_1 = 1, \lambda_2 = -1$。

$\lambda = 1$对应的特征向量即为全体实对称矩阵$\mathbb{S}^2$；$\lambda = -1$对应的特征向量为全体反对称矩阵。那么，若以$E_{11}, E_{22}, E_{12}+E_{21}, E_{12}-E_{21}$为基，则$T$对应的矩阵
$$
T = \begin{pmatrix}
    1 & 0 & 0 & 0\\
    0 & 1 & 0 & 0\\
    0 & 0 & 1 & 0\\
    0 & 0 & 0 & -1\\
\end{pmatrix}
$$

**例7** 证明$\mathbb{S}^2 \oplus \bar{S}^2 = \mathbb{R}, \bar{S}^2 = \mathbf{span}\{E_{12} - E_{21}\}$。

由这两个线性空间的基线性无关可证。也可根据定义证明，即任取$A\in \mathbb{R}^{2\times 2}$，则
$$
A = \frac{1}{2}(A + A^T) + \frac{1}{2}(A - A^T)
$$
又由于$\mathbb{S}^2 \cap \bar{S}^2 = {0}$，因此得证。

**例8** 证明$\mathbb{S}^2$与$\bar{S}^2$均为$\mathcal{T}$的不变子空间。

根据定义，可证$\mathcal{T}A \in V, A\in V$。

这里有一个更一般的结论，若$\mathcal{A}$为$V$上的线性变换，$\lambda$为特征值，令$E_\lambda(\mathcal{A})$为$\lambda$对应的特征向量与零构成的子集，则$E_\lambda(\mathcal{A})$为$V$关于$\mathcal{A}$的不变子空间，称为**特征子空间**。

**证明** $\mathcal{A}(v_1 + v_2) = \mathcal{A}v_1 + \mathcal{A}v_2 = \lambda v_1 + \lambda v_2 = \lambda(v_1 + v_2)$，即$v_1 + v_2\in E_\lambda(\mathcal{A})$；类似地，$kv\in E_\lambda(\mathcal{A})$。又由于$\mathcal{A}v = \lambda v \in E_\lambda(\mathcal{A})$，因此得证。

**定理** 若$\mathcal{A}$为$V$上的线性变换，$\mathcal{A}$可对角化当且仅当$V$可写为特征子空间的直和。