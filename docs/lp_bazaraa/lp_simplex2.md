# 初始解与收敛性

## 初始基可行解

前一章所用的单纯形法是建立在我们知道初始基可行解的基础上的，即我们已知$\mathbf{x}$满足$\mathbf{Ax} = \mathbf{b}, \mathbf{x}\geq \mathbf{0}$，但是在一些情况下这样的初始基可行解是比较难找到的。这一节的目的就是帮我们快速找到初始基可行解。

我们先看这几种约束情况：
- $\mathbf{Ax} \leq \mathbf{b}, \mathbf{x}\geq \mathbf{0}$且$\mathbf{b}\geq \mathbf{0}$。这时，我们通过添加松弛变量可以将约束转化为$\mathbf{Ax} + \mathbf{x}_s = \mathbf{b}, \mathbf{x}\geq \mathbf{0}, \mathbf{x}_s\geq \mathbf{0}$，于是$(\mathbf{x}^T, \mathbf{x}_s^T) = (\mathbf{0}^T, \mathbf{b}^T)$即为一组基可行解，我们从这里开始单纯形法的迭代。
- $\mathbf{Ax} \leq \mathbf{b}, \mathbf{x}\geq \mathbf{0}$，但$\mathbf{b}\not\geq \mathbf{0}$。我们仍然可以通过添加松弛变量可以将约束转化为$\mathbf{Ax} + \mathbf{x}_s = \mathbf{b}, \mathbf{x}\geq \mathbf{0}, \mathbf{x}_s\geq \mathbf{0}$，但这时$\mathbf{b}\not\geq \mathbf{0}$，因此$(\mathbf{x}^T, \mathbf{x}_s^T) = (\mathbf{0}^T, \mathbf{b}^T)$*不可行*。
- 同样，当约束为$\mathbf{Ax} \geq \mathbf{b}, \mathbf{x}\geq \mathbf{0}$且$\mathbf{b}\not\leq \mathbf{0}$时，即使我们将约束两边同时乘以$-1$，但也不能保证约束右边非负。

总的来说，根据第一章的内容，我们知道一个线性规划问题一定能转化为
$$
\begin{array}{lrrrr}
\operatorname{minimize} & \mathbf{c}^T\mathbf{x} \\
\text {subject to} & \mathbf{Ax} & = & \mathbf{b}\\
& \mathbf{x} & \geq & \mathbf{0}
\end{array}
$$
其中$\mathbf{b}\geq 0$，如果$m\times n$阶矩阵$\mathbf{A}$中包含一个子阵是$m$阶单位阵，那么初始基可行解很容易就能得到；但如果没有这样的单位阵，那么初始基可行解很难一下子找到。注意，这里不仅要找一个*基*，也要它是*可行*的！显然不能一个一个尝试，特别是当线性规划问题不可行的时候。

于是，我们引入了**人工变量**来帮助我们解决这些问题。我们可以将原问题转化为
$$
\begin{array}{lrrrrrrr}
\operatorname{minimize} & \mathbf{c}^T\mathbf{x} \\
\text {subject to} & \mathbf{Ax} & + & \mathbf{I}\mathbf{x}_a & = & \mathbf{b}\\
& \mathbf{x} & , & \mathbf{x}_a & \geq & \mathbf{0}
\end{array}
$$
这时，$(\mathbf{x}^T, \mathbf{x}_a^T) = (\mathbf{0}^T, \mathbf{b}^T)$即为一组基可行解。

当然，这里也可以仅引进一个人工变量。我们若找到一组基解（但可能不可行，不要求可行就很容易，只需要高斯消元即可），那么则有
$$
\mathbf{x}_{B}+\mathbf{B}^{-1} \mathbf{N} \mathbf{x}_{N}=\mathbf{B}^{-1} \mathbf{b}
$$
若$\mathbf{B}^{-1} \mathbf{b}\not\geq 0$，则说明该基解不可行。此时，我们添加一个人工变量$x_a$，并令向量$\mathbf{a}_a$满足
$$
({\mathbf{a}_a})_i = \begin{cases}
    0, (\mathbf{B}^{-1} \mathbf{b})_i \geq 0\\
    -1, \text{otherwise}
\end{cases}
$$
接着，改写约束为
$$
\mathbf{x}_{B}+\mathbf{B}^{-1} \mathbf{N} \mathbf{x}_{N} + \mathbf{a}_ax_a=\mathbf{B}^{-1} \mathbf{b}
$$
这时，以$x_a$为进基变量，以第$r$行的基变量为出基变量，其中
$$
r = \operatorname{arg} \min\limits_{i:(\mathbf{B}^{-1} \mathbf{b})_i < 0} \{(\mathbf{B}^{-1} \mathbf{b})_i\}
$$
进行一次换基后，等式右边均会为非负数，我们也就找到了带有一个人工变量的线性规划问题的初始基可行解。

但是，需要注意的是，添加人工变量改变了原问题。只有人工变量最后都为$0$，才有$\mathbf{Ax} = \mathbf{b}$。所以，人工变量只是我们得到初始基可行解的手段，“兔死狗烹，鸟尽弓藏”，发挥完它的作用后，要“迫使”它变为$0$，恢复该线性规划问题最初的面貌。下面将介绍扔掉人工变量的两种方法：两阶段法和大M法。

## 两阶段法

顾名思义，两阶段法包括两个阶段。在第一个阶段中，我们求解线性规划问题
$$
\begin{array}{lrrrrrrr}
\operatorname{minimize} &&& \mathbf{1}^T\mathbf{x}_a \\
\text {subject to} & \mathbf{Ax} & + & \mathbf{I}\mathbf{x}_a & = & \mathbf{b}\\
& \mathbf{x} & , & \mathbf{x}_a & \geq & \mathbf{0}
\end{array}
$$
这一阶段的主要目的有两个，一是“扔掉”人工变量，二是找到替代人工变量的另一基可行解。

如果扔不掉人工变量，意味着上述线性规划问题的最优值不是$0$，于是原问题不可行。假设原问题有可行解$\mathbf{x}$，那么$(\mathbf{x}^T, \mathbf{x}_a^T) = (\mathbf{x}^T, \mathbf{0}^T)$对于添加人工变量后的线性规划问题也是可行的，此时目标函数值为$0$。

如果原问题可行，那么上述线性规划问题取最优值时，所有的人工变量应为$0$。在没有退化发生的情况下，人工变量应该都不是基变量，因此我们就找到了原问题的一组基可行解。在第二个阶段中，我们只要对原问题以这组基可行解为始，进行单纯形法的迭代即可。

因此，使用两阶段法进行求解有下列几种情形：

### 情形1

第一阶段解得$\mathbf{x}_a\neq {0}$，此时原线性规划问题无解。

### 情形2

第一阶段解得$\mathbf{x}_a = {0}$，且迭代终止时$\mathbf{x}_a$为非基变量。此时，第一阶段迭代完成后，单纯形表应为

$$
\begin{array}{c|c|ccc|c|c}
& x_0 & \mathbf{x}_B & \mathbf{x}_N & \mathbf{x}_a &\text{RHS} & \\\hline
x_0 & 1 & \mathbf{0} & \mathbf{0} & -\mathbf{1} & 0 & \text{Row 0} &\\
\mathbf{x}_B & \mathbf{0} & \mathbf{I} & \mathbf{B}^{-1} \mathbf{N} & \mathbf{B}^{-1} & \mathbf{B}^{-1}  \mathbf{b} & \text{Rows 1 through m} &
\end{array}
$$

注意，一般情况下的单纯形表如下，只是这里$\mathbf{c}_{B} = \mathbf{c}_{N} = \mathbf{0}, \mathbf{c}_a = \mathbf{1}, \mathbf{N}_a = \mathbf{I}$。通过第一阶段，我们找到了一组基可行解。丢弃人工变量，从当前这组基开始迭代。其中，$\sigma_j$和$z_0$可以通过$\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{N}- \mathbf{c}_{N}$和$\mathbf{c}_{B} \mathbf{B}^{-1}  \mathbf{b}$计算得到（通过初等行变换将$\mathbf{x}_B$的系数消为$0$即可）。

$$
\begin{array}{c|c|cc|c|c}
& z & \mathbf{x}_B & \mathbf{x}_N & \text{RHS} & \\\hline
z & 1 & \mathbf{0} & \mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{N}- \mathbf{c}_{N} & \mathbf{c}_{B} \mathbf{B}^{-1}  \mathbf{b} & \text{Row 0} &\\
\mathbf{x}_B & \mathbf{0} & \mathbf{I} & \mathbf{B}^{-1} \mathbf{N} & \mathbf{B}^{-1}  \mathbf{b} & \text{Rows 1 through m} &
\end{array}
$$

### 情形3

第一阶段解得$\mathbf{x}_a = {0}$，且迭代终止时$\mathbf{x}_a$中存在值为$0$的基变量。此时的单纯形表应为
$$
\begin{array}{c|c|cccc|c|c}
& x_0 & \mathbf{x}_B & \mathbf{x}_N & \mathbf{x}_{a_N} & \mathbf{x}_{a_B} &\text{RHS} & \\\hline
x_0 & 1 & \mathbf{0} & \mathbf{0} & -\mathbf{1} & \mathbf{0} & 0 & \text{Row 0} &\\
\mathbf{x}_B & \mathbf{0} & \mathbf{I} & \mathbf{R}_1 & \mathbf{R}_3 & \mathbf{0} & \mathbf{\bar{b}} & \text{Rows 1 through k} &\\
\mathbf{x}_{a_B} & \mathbf{0} & \mathbf{0} & \mathbf{R}_2 & \mathbf{R}_4 & \mathbf{I} & \mathbf{0} & \text{Rows k+1 through m} &
\end{array}
$$

我们处理这种情况的方法之一就是把人工变量“赶出”基。例如如果我们想赶出$({x}_{a_{B}})_i$，只要矩阵$\mathbf{R}_2$的第$i$行不全为$0$即可，从这一行里选择不为$0$的一列，该列对应的变量作为进基变量，这也满足前面所说进出基需要的三个要求。
- 进出基变换后，所有基变量还是一组基（即系数矩阵满秩）；
- 求得的基解可行（即$\mathbf{x}_B \geq 0$）；
- 目标函数值减少或不变化（对于最小化问题）。

但若$\mathbf{R}_2$的某一行全为$0$，这就无法选择进基变量使得其满足第一个要求了。这是由于这一行的约束是“多余”的（即与其他行间线性相关），因此可以直接删去这条约束、删去这一行。

依次一步一步，最终我们可以将所有的人工变量全部变为非基变量。接着，就可以开启我们的第二阶段了。

另一种方法是直接进入第二阶段，我们可以根据Dantzig规则选取进基变量，但是在选取出基变量时，不能使用之前的规则进行出基。因为
$$
\mathbf{x}_{B} = \begin{pmatrix}
    \bar{b_1}\\ \vdots \\ \bar{b_m}
\end{pmatrix} - 
\begin{pmatrix}
    {y_{1k}}\\ \vdots \\ {y_{mk}}
\end{pmatrix}x_k \geq 0
$$
若$\mathbf{R}_2$中进基变量对应的这列有负值，那么选择其他变量出基后，本该为$0$的人工基变量就会增加，就不再是$0$。因此，我们应该以这些负值所在行对应的人工基变量作为出基变量，一旦出基，就“踢走”这个人工变量。这是对原来迭代过程的一点小小改变。

## 大M法

在大M法中，我们引入了一个很大的数$M$（至于$M$应该有多大，$\mathbf{c}$和$\mathbf{A}$都应该关注，特别是当$\mathbf{A}$很病态时），求解线性规划问题$P(M)$：
$$
\begin{array}{lrrrrrrr}
\operatorname{minimize} &\mathbf{c}^T\mathbf{x}&+& M\mathbf{1}^T\mathbf{x}_a \\
\text {subject to} & \mathbf{Ax} & + & \mathbf{I}\mathbf{x}_a & = & \mathbf{b}\\
& \mathbf{x} & , & \mathbf{x}_a & \geq & \mathbf{0}
\end{array}
$$

这个问题显然可以找到一个基可行解$\mathbf{x} = \mathbf{0}, \mathbf{x}_a = \mathbf{b}$。另外，如果原问题有可行解，所有的人工变量最后都会被强制变为$0$。这是因为有一个很大的系数$M$在各人工变量前面，这就相当于让$\mathbf{1}^T\mathbf{x}_a$最小是求解这个问题的首要任务，这个任务的权重很大（有$M$这么大），这个问题完成不好，最小化$\mathbf{c}^T\mathbf{x}$完成得好也无济于事。那这个首要任务怎么才算完成得好呢，就是尽可能使$\mathbf{x}_a = \mathbf{0}$，让这个问题回归本来的面貌。因此，大M法本质上就是两阶段法，只不过现在目标函数变为了两个阶段目标函数的加权和，也就是：
$$
z_{\text{big-M}} = Mz_{\text{phase I}} + z_{\text{phase II}}\\
z_{\text{phase I}} = \mathbf{1}^T\mathbf{x}_a\\
z_{\text{phase II}} = \mathbf{c}^T\mathbf{x}
$$

下面，我们关心的还是原问题什么时候有最优解、什么时候最优解无界和什么时候无可行解。根据两阶段法，第一阶段结束后可以根据人工变量是否为$0$判断原问题是否可行，这通过$z_{\text{big-M}} = Mz_{\text{phase I}} + z_{\text{phase II}}$反映在大M法上也就是，当所有$\sigma_j$中$M$的系数都非正时（这对应着两阶段法中第一阶段里所有$\sigma_j \leq 0$），判断人工变量这个时候是不是都是$0$。如果人工变量都是$0$，那么我们就通过大M法找到了一个基可行解；否则，说明该问题不可行，这个时候迭代应该停止。当所有人工变量都是$0$时，一方面我们可以像两阶段法那样，扔掉人工变量（可以参考前面两阶段法中介绍的如何扔掉基变量中等于$0$的人工变量）；另一方面，在大M法中，我们也可以选择继续正常迭代，直到最优解存在或者无界的情况。

## 退化

前面说到，如果在极点处有多于$n$个的线性无关的“紧”约束，该极点对应的基可行解$\mathbf{x}_B$中就有值为$0$，此时就称发生了退化。退化可能会导致一些不同寻常的结果，包括：
- 当检验数$\mathbf{c}_{B}  \mathbf{B}^{-1} \mathbf{N} - \mathbf{c}_{N}$中有$0$时，不一定意味着该问题有无数个最优解（这在没有退化时是正确的），有可能最优解对应的极点是退化的，所以导致$\mathbf{c}_{B}  \mathbf {B}^{- 1 }  \mathbf {N}- \mathbf {c}_{N}\not< 0$，此时这个检验数$\sigma_j = 0$的非基变量也是一点儿也不能增大的。
- 单纯形法的迭代可能“进循环了”。退化可能会导致这样一种情况，单纯形法中的“换基”操作换来换去还是一个极点，于是一不小心，基换了一轮又回到了原来的某个基。这往往是由于在每一轮迭代中，因为退化而导致的有多个候选的出基变量可供选择。这个时候，就需要通过一些方法避免循环。循环在没有退化时不会发生，是因为没有退化时，目标函数的值总在不断减少的，也就意味着每次迭代的基都不会跟之前迭代过的基重复。

下面将介绍两种避免循环的迭代规则，并说明通过这些规则，也会实现像没有退化时一样的**没有一个基会被重复迭代到**。同时，由于基的数量是有限的，因此在这些规则下，迭代总能终止。但是，需要注意的是，即使我们有办法通过这些规则来避免循环，由于退化的存在，单纯形法仍然有可能在某一个极点处迭代很多次，导致单纯形法无法在多项式的时间复杂度内完成计算。同时，由于这些避免循环的规则也需要较多的时间进行，所以商业求解器中很少应用这些规则。

### 字典规则

在该规则下，定义
$$
I_{0}=\left\{r: \frac{\bar{b}_{r}}{y_{r k}}=\min_{1 \leq i \leq m}\left\{\frac{\bar{b}_{i}}{y_{i k}}: y_{i k}>0\right\}\right\}\\
I_{j}=\left\{r: \frac{y_{r j}}{y_{r k}}=\min_{i \in I_{j-1}}\left\{\frac{y_{i j}}{y_{i k}}\right\}\right\}
$$
其中，$k$是选中的进基变量的序号，$y_{ij}$是单纯形表中每次迭代的系数。如果集合$I_0$里面只有一个元素，那么这个元素就是我们这一轮应该选择的出基变量的序号；否则就考察集合$I_1$，以此类推。到了集合$I_m$肯定只有一个元素，否则，$\mathbf{B}^{-1}$就会有至少两行是成比例的，也就是线性相关的，这和$\mathbf{B}^{-1}$满秩构成了矛盾。从计算时间上来看，这样的规则可能导致每一次迭代花费较长的时间。

为了证明这种规则能够避免循环，我们首先定义向量$\mathbf{x}\succ 0$为
- $\mathbf{x}\neq \mathbf{0}$；
- $\mathbf{x}$的第$1$个非零元素为正。

在字典规则中，假设单纯形法从如下所示的情形开始迭代，$\mathbf{x}_0$作为初始基变量。

$$
\begin{array}{c|c|cc|c|c}
& z & \mathbf{x}_0 & \mathbf{x}_1 & \text{RHS} & \\\hline
z & 1 & \mathbf{0} & \mathbf{c}_{B}\mathbf{A}- \mathbf{c}_{1} & \mathbf{c}_{B}\mathbf{b} & \text{Row 0} &\\
\mathbf{x}_0 & \mathbf{0} & \mathbf{I} &  \mathbf{A} & \mathbf{b} & \text{Rows 1 through m} &
\end{array}
$$

首先可以证明根据字典规则，每一轮迭代后矩阵$(\mathbf{B}^{-1}\mathbf{b}, \mathbf{B}^{-1})$（也就是上表中**RHS**一列与$\mathbf{x}_0$对应的$m$列向量组成的矩阵）的每一行形成的向量$\succ 0$。初始时该结论显然成立。每一次迭代时，假设选择了第$r$行第$k$列的元素作为主元，则换基后该矩阵第$i$行为
$$
\left(\frac{\bar{b}_{r}}{y_{r k}}, \frac{y_{r 1}}{y_{r k}}, \ldots, \frac{y_{r m}}{y_{r k}}\right), \text{if }i = r\\
\left(\bar{b}_{i}, y_{i 1}, \ldots, y_{i m}\right)-\frac{y_{i k}}{y_{r k}}\left(\bar{b}_{r}, y_{r 1}, \ldots, y_{r m}\right), \text{if } i \neq r
$$
当$i = r$时，显然向量依然满足$\succ 0$；当$i \neq r$时，需要依次分情况讨论。$y_{ik}\leq 0$时，向量满足结论；否则就需要讨论$\bar{b}_{i}$与$\frac{y_{i k}}{y_{r k}}\bar{b}_{r}$的大小，而这有两种情况（大于或等于，只要根据$i$是否属于$I_0$而判断）。若大于，则向量满足结论；否则，就要继续讨论$y_{i1}$与$\frac{y_{i k}}{y_{r k}}y_{r1}$的大小，而这又是根据$i$是否属于$I_1$而判断的。……以此类推，结论成立。

迭代中的单纯形表如下所示。
$$
\begin{array}{c|c|cc|c|c}
& z & \mathbf{x}_0 & \mathbf{x}_1 & \text{RHS} & \\\hline
z & 1 & \mathbf{c}_{B}\mathbf{B}^{-1} - \mathbf{c}_0 & \mathbf{c}_{B}\mathbf{B}^{-1}\mathbf{A}- \mathbf{c}_{1} & \mathbf{c}_{B}\mathbf{B}^{-1}\mathbf{b} & \text{Row 0} &\\
\mathbf{x}_0 & \mathbf{0} & \mathbf{B}^{-1} &  \mathbf{B}^{-1}\mathbf{A} & \mathbf{B}^{-1}\mathbf{b} & \text{Rows 1 through m} &
\end{array}
$$

当一次迭代基从$\mathbf{B}$换至$\hat{\mathbf{B}}$时，有
$$
\left(\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{b}, \mathbf{c}_{B} \mathbf{B}^{-1}\right)-\left(\mathbf{c}_{\hat{B}} \hat{\mathbf{B}}^{-1} \mathbf{b}, \mathbf{c}_{\hat{B}} \hat{\mathbf{B}}^{-1}\right)=\frac{\sigma _{k}}{y_{r k}}\left(\bar{b}_{r}, y_{r 1}, y_{r 2}, \ldots, y_{r m}\right)
$$
由于$x_k$是出基变量，因此$\sigma_k > 0, y_{rk} > 0$，因此
$$
\left(\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{b}, \mathbf{c}_{B} \mathbf{B}^{-1}\right)-\left(\mathbf{c}_{\hat{B}} \hat{\mathbf{B}}^{-1} \mathbf{b}, \mathbf{c}_{\hat{B}} \hat{\mathbf{B}}^{-1}\right) \succ 0
$$

假设碰到循环了，那么经过一次循环，我们会得到
$$
\left(\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{b}, \mathbf{c}_{B} \mathbf{B}^{-1}\right) - \left(\mathbf{c}_{B} \mathbf{B}^{-1} \mathbf{b}, \mathbf{c}_{B} \mathbf{B}^{-1}\right) = \mathbf{0}\succ 0
$$
这就是矛盾所在。

### Bland规则

Bland规则同时作用于进基和出基变量的选择，它听起来比Danzig规则简单。首先，为所有的变量设定一个优先级（不妨取$x_1, x_2, \cdots, x_n$）；接着，对于进基变量的选择，在所有检验数$\mathbf{c}_{B}  \mathbf{B}^{-1} \mathbf{N} - \mathbf{c}_{N} > 0$的非基变量里，选择优先级最高的进基；最后，对于出基变量，当可供选择的出基变量不止一个时，也选择优先级最高的出基。从计算时间上来看，这样的规则对于单次迭代可能影响不大，但是相较于Danzig规则，它的迭代次数可能较多。

我们使用反证法证明在Bland规则下不会发生循环。假设发生了循环，基经历了$B_1\rightarrow B_2 \rightarrow \cdots \rightarrow B_l\rightarrow B_1$（不妨设循环从第$1$轮迭代开始）。在这个过程中肯定有变量出基完又在后续迭代中重新入基，在所有的这些变量中，我们考察**优先级最低**的变量$x_t$（因此，下文中的变量$x_s$和$x_r$的优先级均比$x_t$高），假设$x_t$在第$k$轮迭代出基了，接着又重新在第$f$轮迭代入基，它出基的同时入基的变量为$x_s$。

在第$k$轮迭代时，有方程组（1）
$$
\mathbf{x}_{B_k} + \mathbf{N}_k\mathbf{x}_{N_k} = \mathbf{b}_k\\
z = z_0 + \mathbf{c}^T_k \mathbf{x}_{N_k}
$$
由于选择$x_s$变量入基，因此
$$
(\mathbf{c}_k)_s < 0, ({\mathbf{a}_k}_s)_t > 0, (\mathbf{b}_k)_t = 0
$$

当循环发生时，目标函数值一直没有改变（因为迭代不可能使目标函数值变大，而换基换了一圈又换回了起点，相同的基对应的目标函数值是一样的）。在第$f$轮迭代时，有方程组（2）
$$
\mathbf{x}_{B_f} + \mathbf{N}_f\mathbf{x}_{N_f} = \mathbf{b}_f\\
z = z_0 + \mathbf{c}^T_f \mathbf{x}_{N_f}
$$
此时，由于$x_t$重新入基，因此
$$
(\mathbf{c}_f)_t < 0
$$
$x_s$优先级比$x_t$高，但根据Bland规则$x_s$没有入基，说明
$$
(\mathbf{c}_f)_s \geq 0
$$
（$x_s$为基变量时，$(\mathbf{c}_f)_s = 0$）。

因为单纯形法的迭代仅仅只是各方程之间的加加减减（代数意义上也就是系数矩阵的初等行变换），满足方程组（1）的解也一定满足方程组（2）（即使可能会不满足约束$\mathbf{x}\geq 0$），于是我们将方程组（1）的一个解
$$
(\mathbf{x})_s = q, (\mathbf{x})_j = 0(j\in N_k \wedge j\neq s), \mathbf{x}_{B_k} = \mathbf{b}_k - {\mathbf{a}_k}_s q, \\
z = z_0 + (\mathbf{c}_k)_s q
$$
代入方程组（2），得到
$$
z = z_0 + (\mathbf{c}_k)_sq = z_0 + (\mathbf{c}_f)_s q + \sum_{x_j\in B_k} (\mathbf{c}_f)_j(\mathbf{b}_k - {\mathbf{a}_k}_s q)_j
$$

因此
$$
\left[(\mathbf{c}_k)_s - (\mathbf{c}_f)_s + \sum_{x_j \in B_k}(\mathbf{c}_f)_j({\mathbf{a}_k}_s)_j\right]\cdot q \equiv \sum_{x_j \in B_k}(\mathbf{c}_f)_j({\mathbf{b}_k})_j\\
\Rightarrow (\mathbf{c}_k)_s - (\mathbf{c}_f)_s + \sum_{x_j \in B_k}(\mathbf{c}_f)_j({\mathbf{a}_k}_s)_j = \sum_{x_j \in B_k}(\mathbf{c}_f)_j({\mathbf{b}_k})_j = 0
$$
由于$(\mathbf{c}_k)_s < 0, (\mathbf{c}_f)_s \geq 0$，因此
$$
\sum_{x_j \in B_k}(\mathbf{c}_f)_j({\mathbf{a}_k}_s)_j > 0
$$
因此，存在$x_r\in B_k$满足
$$
(\mathbf{c}_f)_r({\mathbf{a}_k}_s)_r > 0
$$
可以看到，$(\mathbf{c}_f)_r \neq 0\Rightarrow x_r\not\in B_f$，因此变量$x_r$也是在该循环里出基完又入基，因此它的优先级也比$x_t$高。但是，它在第$f$轮迭代中没有入基，说明
$$
(\mathbf{c}_f)_r \geq 0\Rightarrow ({\mathbf{a}_k}_s)_r > 0, (\mathbf{c}_f)_r > 0
$$
综上所述，$x_r$在第$k$轮迭代时也可以被选作出基变量（$({\mathbf{a}_k}_s)_r > 0, (\mathbf{b}_k)_r = 0$），并且它的优先级比$x_t$高，根据Bland规则，它应该被选作出基变量，这与前面$x_t$被选作出基变量矛盾。所以，假设不成立，Bland规则下不会发生循环。