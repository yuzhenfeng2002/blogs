# 离散仿真优化——局部收敛算法简介

## COMPASS算法

COMPASS算法的特点包括：
1. 拥有自适应邻域结构；
2. 能够求解完全约束问题、部分约束问题、无约束问题；
3. 收敛至局部最优解的集合。

### 完全约束下的COMPASS算法

考虑以下的离散仿真优化问题：
$$
\min\limits_{x\in \Theta} E_{\psi}[G(\mathbf{x}, \psi)]
$$
其中，$\Theta = \Phi \cap \mathcal{Z}^d$，其中$\Phi$是闭且有界的，$\mathcal{Z}^d$是$d$维整数集合；$\psi$是仿真的随机输入，用$G_i(\mathbf{x})$表示对$G(\mathbf{x}, \psi)$的第$i$次观测。假设$G(\mathbf{x}, \psi)$是可测且可积的，记$g(\mathbf{x}) = E_{\psi}[G(\mathbf{x}, \psi)]$。

> 假设 1
> 
> 对于每一$\mathbf{x}\in\Theta$，有
> $$
  P\{\lim\limits_{r\rightarrow\infty}\frac{1}{r}\sum\limits_{i=1}^rG_i(\mathbf{x}) = g(\mathbf{x})\} = 1
> $$

#### 算法

```
-------------------- INITIALIZE --------------------
k = 0
find x[0] from Theta
V[k] = {x[0]}
take SAR(k, x) observations and calculate the sample mean AG[x] for every x in V[k]
C[k] = Theta
-------------------- ITERATIONS --------------------
while stopping condition == False:
    k = k + 1
    sample m elements y[1], ..., y[m] from C[k - 1]
    V[k] = V[k - 1] U {y[1], ..., y[m]}
    take SAR(k, x) observations and calculate the sample mean AG[x] for every x in V[k]
    xmin = argmin(AG)
    C[k] = MPA(xmin, V[k])
```

其中，`V[k]`指经过$k$次迭代已经访问过的解集，用$\mathcal{V}_k$表示；`SAR(k, x)`指第$k$次迭代时为解$\mathbf{x}$分配的仿真次数，用$a_k(\mathbf{x})$表示；`AG[x]`指$k$次迭代中在该解处的平均值，用$\bar{G}_k(\mathbf{x})$表示；`C[k]`指第$k$次迭代后构造的最有前景的区域（most promising area），用$\mathcal{C}_k$表示；`xmin`指当前访问过的解集中具有最小平均值的解，用$\hat{\mathbf{x}}_k^\star$表示；`MPA(xmin, V[k])`指根据当前结果构造最有前景的区域的规则，在该算法中$\mathcal{C}_{k}=\left\{\mathbf{x}: \mathbf{x} \in \Theta \text{ and }\left\|\mathbf{x}-\hat{\mathbf{x}}_{k}^{\star}\right\| \leq \|\mathbf{x}-\mathbf{y}\| \forall \mathbf{y} \in \mathcal{V}_{k}\text{ and }\mathbf{y} \neq \hat{\mathbf{x}}_{k}^{\star}\right\}$。

#### 收敛性

> 假设 2
> 
> $a_k(\mathbf{x})$满足：
> 1. $a_k(\mathbf{x}) \geq 1 \text{ if } \mathbf{x}\in \mathcal{V}_k-\mathcal{V}_{k-1}$；
> 2. $\lim\limits_{k \rightarrow \infty} N_{k}(\mathbf{x})=+\infty$，对于所有$\mathbf{x} \in \bigcup_{k=0}^{\infty} \mathcal{V}_{k}$，其中$N_{k}(\mathbf{x})$为$k$次迭代共计在解$\mathbf{x}$处进行仿真的次数。

> 定理 1
> 
> 在满足假设1与2的情况下，序列$\{\hat{\mathbf{x}}_k^\star\}_{k=0}^{\infty}$满足$\operatorname{Pr}\{\hat{\mathbf{x}}_{k}^{\star} \notin \mathcal{M}\text{ i.o.}\}=0$，其中$\mathcal{M} = \{\mathbf{x}: \mathbf{x}\in \Theta\text{ and either }\mathcal{N}(\mathbf{x})=\varnothing\text{ or }g(\mathbf{x}) \leq g(\mathbf{y})\text{ for all }\mathbf{y} \in \mathcal{N}(\mathbf{x})\}, \mathcal{N}(\mathbf{x}) = \{\mathbf{y}: \mathbf{y}\in \Theta\text{ and }\|\mathbf{x} - \mathbf{y}\| = 1\}$。

> 补充 **i.o.**
> 
> Let $\left\{A_{n}\right\}_{n=1}^{\infty}$ be an infinite sequence of events. We say that events in the sequence occur "infinitely often" if $A_{n}$ holds true for an infinite number of indices $n \in\{1,2,3, \ldots\} .$ We say that events in the sequence occur "finitely often" if they do not occur infinitely often, that is, if $A_{n}$ holds true for at most finitely many indices $n \in\{1,2,3, \ldots\} .$ We formally define the set of all outcomes $\omega$ for which an infinite number of the events is true, and for which at most finitely many events are true, as follows:
> $$
  \left\{A_{n} \text { i.o. }\right\} =\cap_{m=1}^{\infty} \cup_{n=m}^{\infty} A_{n} \\
  \left\{A_{n} \text { f.o. }\right\} = \{A_{n} \text { i.o. }\}^{c} = \cup_{m=1}^{\infty} \cap_{n=m}^{\infty} A_{n}^{c}
> $$
> where the final equality holds by DeMorgan's law.
> 
> The right-hand-side of (1) can be read as follows: For all positive integers $m$, there exists a positive integer $n \geq m$ such that $A_{n}$ is true. Some thought will convince you that this holds if and only if an infinite number of the events are true. The right-hand-side of (2) can be read as follows: There is a positive integer $m$ such that $A_{n}$ is false for all integers $n \geq m$. Some thought will convince you that this holds if and only if an at most finite number of the events occur.

#### 证明思路

首先对要证明的定理进行等价转换：
$$
\operatorname{Pr}\{\hat{\mathbf{x}}_{k}^{\star} \notin \mathcal{M}\text{ i.o.}\} = 0\\
\Leftrightarrow \sum_{A\subset \Theta} \operatorname{Pr}\{\hat{\mathbf{x}}_{k}^{\star} \notin \mathcal{M}\text{ i.o.} | V_{\infty} = A\}\operatorname{Pr}\{V_{\infty} = A\} =0\\
\Leftrightarrow \operatorname{Pr}\{\hat{\mathbf{x}}_{k}^{\star} \notin \mathcal{M}\text{ i.o.} | V_{\infty} = A\} = 0 \text{ for all nonempty set } A \subset \Theta \text{ such that } \operatorname{Pr}\{V_{\infty} = A\} > 0
$$

若$\hat{\mathbf{x}}_{k}^{\star} \notin \mathcal{M}\text{ i.o.}$，则存在$\mathbf{x}\in A \text{ and } \mathbf{x}\notin \mathcal{M}$ 满足$\hat{\mathbf{x}}_{k}^{\star} = \mathbf{x}\text{ i.o.}$。此时，一定存在$\mathbf{y} \in \mathcal{N}(\mathbf{x})$满足$g(\mathbf{x}) > g(\mathbf{y})$，则$\mathbf{y}\in \mathcal{C}_{k}$，因此，
$$
\operatorname{Pr}\{\mathbf{y}\in \mathcal{V}_{k+1} | \hat{\mathbf{x}}_{k}^{\star} = \mathbf{x}, \mathbf{y}\notin \mathcal{V}_{k}\} \geq \frac{1}{|\Theta|} > 0\\
\Rightarrow \operatorname{Pr}\{\mathbf{y}\in A | \mathcal{V}_{\infty} = A, \hat{\mathbf{x}}_{k}^{\star} = \mathbf{x}\text{ i.o.}\} = 1
$$
这意味着$g(\hat{\mathbf{x}}_{k}^{\star}) \neq \min_{\mathbf{z}\in A} g(\mathbf{z}) \text{ i.o.}$，下证$P\{\lim\limits_{k\rightarrow \infty} g(\hat{\mathbf{x}}_{k}^{\star}) = \min_{\mathbf{z}\in \mathcal{V}_{\infty}} g(\mathbf{z}) | \mathcal{V}_{\infty} = A\} = 1$，由此可以导出矛盾。

根据三角不等式，
$$
\operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\min _{\mathbf{x} \in \mathcal{V}_{\infty}} g(\mathbf{x})\right| \geq \epsilon \text { i.o.} \mid \mathcal{V}_{\infty}=A\right\}\\
\leq \operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\bar{G}_{k}\left(\hat{\mathbf{x}}_{k}^{\star}\right)\right| \geq \frac{\epsilon}{2} \text{ i.o.} \mid \mathcal{V}_{\infty}=A\right\} + \operatorname{Pr}\left\{\left|\bar{G}_{k}\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\min _{\mathbf{x} \in \mathcal{V}_{\infty}} g(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text{ i.o.} \mid \mathcal{V}_{\infty}=A\right\}\\
\leq \operatorname{Pr}\left\{\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text { i.o. for some } \mathbf{x} \in \mathcal{V}_{\infty} \mid \mathcal{V}_{\infty}=A\right\} + \operatorname{Pr}\left\{\left|\min _{\mathbf{x} \in \mathcal{V}_{\infty}} \bar{G}_{k}(\mathbf{x})-\min _{\mathbf{x} \in \mathcal{V}_{\infty}} g(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text { i.o.} \mid \mathcal{V}_{\infty}=A\right\}
$$

根据引理，
$$
\leq \operatorname{Pr}\left\{\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text { i.o. for some } \mathbf{x} \in \mathcal{V}_{\infty} \mid \mathcal{V}_{\infty}=A\right\} + \operatorname{Pr}\left\{\max _{\mathbf{x} \in \mathcal{V}_{\infty}}\left|\bar{G}_{k}(\mathbf{x})-g(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text { i.o.} \mid \mathcal{V}_{\infty}=A\right\}\\
\leq 2 \operatorname{Pr}\left\{\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\epsilon}{2}\text{ i.o. for some }\mathbf{x} \in \mathcal{V}_{\infty} \mid \mathcal{V}_{\infty}=A\right\}
$$

根据Bonferroni不等式，
$$
\leq 2 \sum_{\mathbf{x} \in A} \operatorname{Pr}\left\{\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text{ i.o.} \mid \mathcal{V}_{\infty}=A\right\}\\
=2 \sum_{\mathbf{x} \in A} \operatorname{Pr}\left\{\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text{ i.o.} \mid \mathbf{x} \in \mathcal{V}_{\infty}\right\}
$$

由于
$$
\begin{cases}
  P\{\lim\limits_{r\rightarrow\infty}\frac{1}{r}\sum\limits_{i=1}^rG_i(\mathbf{x}) = g(\mathbf{x})\} = 1\\
  \lim\limits_{k \rightarrow \infty} N_{k}(\mathbf{x})=+\infty\\
  \bar{G}_k(\mathbf{x}) = \frac{1}{N_{k}}\sum\limits_{i=1}^{N_{k}}G_i(\mathbf{x})
\end{cases}
\Rightarrow \operatorname{Pr}\left\{\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\epsilon}{2} \text{ i.o.} \mid \mathbf{x} \in \mathcal{V}_{\infty}\right\} = 0
$$
因此得证。

### MPA的构造与抽样

在最有前景的区域中抽样的方法有Mix-D算法（Pichitlamken and Nelson 2003）与Revised Mix-D算法（Hong and Nelson 2006）。这两种算法用于在凸且紧的集合中抽取整数点。

由于
$$
\mathcal{C}_{k}=\left\{\mathbf{x}: \mathbf{x} \in \Theta \text{ and }\left\|\mathbf{x}-\hat{\mathbf{x}}_{k}^{\star}\right\| \leq \|\mathbf{x}-\mathbf{y}\| \forall \mathbf{y} \in \mathcal{V}_{k}\text{ and }\mathbf{y} \neq \hat{\mathbf{x}}_{k}^{\star}\right\}\\
= \left\{\mathbf{x}: \mathbf{x} \in \Theta \text{ and }(\hat{\mathbf{x}}_{k}^{\star} - \mathbf{y})(\mathbf{x} - \frac{\hat{\mathbf{x}}_{k}^{\star} + \mathbf{y}}{2})\geq 0 \forall \mathbf{y} \in \mathcal{V}_{k}\text{ and }\mathbf{y} \neq \hat{\mathbf{x}}_{k}^{\star}\right\}
$$
则若$\Theta$是凸且紧的集合，则$\mathcal{C}_{k}$也是凸且紧的集合。

两种算法均从一个可行解$\mathbf{x}_0$开始，首先令$\mathbf{X}_0 = \mathbf{x}_0$。
- Mix-D算法首先构造了一个紧紧覆盖凸集的超矩形，在第$t$次迭代时，从超矩形内的整数点里抽取$\mathbf{x}$，在$\mathbf{X}_{t-1}$与$\mathbf{x}$的连线上抽取在集合内的整数点作为$\mathbf{X}_{t}$。该算法的缺点在于当凸集是“歪斜”的时候，$\mathbf{x}_0$与$\mathbf{x}$的连线上满足要求的点较少。
- Revised Mix-D算法中，在第$t$次迭代时，过$\mathbf{X}_{t-1}$作平行于任选一坐标轴的平行线，在这条平行线上抽取符合条件的$\mathbf{X}_{t}$。

> 定理 3
> 
> 记$\mathcal{H}$为$\Theta$通过Mix-D（RMD）从$\mathbf{x}_0$可达的整数点的集合，$\{\mathbf{X}_{t}:t \geq 0\}$的渐进分布为$\mathcal{H}$上的均匀分布。

## LCRS算法框架

LCRS是一个框架，其可以在较温和的条件和假设的情况下依概率1收敛至局部最优解的集合，它包括两部分：
- 抽样（sampling）：从解集中选择可行解；
- 估计（estimation）：对解进行仿真，估计目标函数的值。

### 算法框架

#### 伪代码

```
-------------------- INITIALIZE --------------------
k = 0
find x[0] from Theta
s[k] = {x[0]}
S[k] = s[k]
-------------------- ITERATIONS --------------------
while stopping condition == False:
    # Sampling Scheme
    k = k + 1
    determine C[k] and F[k] according to the sampling scheme
    sample m[k] solutions from C[k] using F[k]
    s[k] = removeDuplicates(m[k])
    S[k] = S[k - 1] U s[k]
    # Estimation Scheme
    determine e[k] which is a subset of S[k]
    take SAR(k, x) observations and calculate the sample mean AG[x] for every x in e[k]
    xmin = argmin(AG)
```

其中，`s[k]`指经过第$k$次抽取的解集，用$\mathcal{S}_k$表示；`S[k]`指经过$k$次迭代一共抽取的解集，用$\mathcal{S}(k)$表示；`C[k]`指第$k$次迭代后构造的最有前景的区域（most promising area），用$\mathcal{C}_k$表示；`F[k]`指第$k$次迭代时在最有前景的区域上用于抽样的分布，用$F_k$表示；`e[k]`指经过第$k$次进行评估的解集，用$\mathcal{E}_k$表示；`SAR(k, x)`指第$k$次迭代时为解$\mathbf{x}$分配的仿真次数，用$a_k(\mathbf{x})$表示；`AG[x]`指$k$次迭代中在该解处的平均值，用$\bar{G}_k(\mathbf{x})$表示；`xmin`指当前访问过的解集中具有最小平均值的解，用$\hat{\mathbf{x}}_k^\star$表示。

#### 条件

算法需要满足的条件包括：

> 条件 1
> 
> - 存在确定的有限集序列$\mathcal{B}_k$，满足$\mathcal{C}_k \subset \mathcal{B}_k, \mathcal{B}_1\subset\mathcal{B}_2\subset\cdots$且$\Theta\subset\cup_{k=1}^{\infty}\mathcal{B}_k$；
> - 抽样分布$F_k$满足$\operatorname{Pr}\{\mathbf{x}\in \mathcal{S}_k\}\geq \epsilon, \forall \mathbf{x} \in \mathcal{N}(\hat{\mathbf{x}}_k^\star) - \mathcal{E}(k - 1)$，其中$\mathcal{E}(k) = \subset\cup_{i=0}^{k}\mathcal{E}_k$。

对于完全约束下的仿真优化问题，条件1.1只需令$\mathcal{B}_k = \Theta$即可。

> 条件 2
>
> - $\mathcal{E}_{k}\subset \mathcal{S}(k)$；
> - $\mathcal{E}_{k}$包括$\mathbf{x}_{0}, \mathcal{N}\left(\widehat{\mathbf{x}}_{k-1}^{*}\right) \cap \mathcal{E}(k-1), \mathcal{S}_{k}$；
> - $a_{k}(\mathbf{x})$满足$\min _{\mathbf{x} \in \mathcal{E}_{k}} N_{k}(\mathbf{x}) \geq 1, \forall k=1,2, \cdots$且$\min _{\mathbf{x} \in \mathcal{E}_{k}} N_{k}(\mathbf{x}) \rightarrow \infty$ w.p. 1 as $k \rightarrow \infty$，其中$N_{k}(\mathbf{x})$为$k$次迭代共计在解$\mathbf{x}$处进行仿真的次数；
> - $|\mathcal{E}(\infty)|<\infty$ w.p. 1。

对于完全约束下的仿真优化问题，条件2.4是一定满足的。

### LCRS算法的局部收敛性

> 假设 1
> 
> 存在$\delta_0 > 0$，使得集合$\mathcal{L} = \{\mathbf{x}\in \Theta: g(\mathbf{x}) \leq g(\mathbf{x}_0) + \delta_0\}$是有限集。

同样，对于完全约束下的仿真优化问题，假设1是一定满足的。

> 假设 2
> 
> 对于每一$\mathbf{x}\in\Theta$，有
> $$
  \operatorname{Pr}\{\lim\limits_{r\rightarrow\infty}\frac{1}{r}\sum\limits_{i=1}^rG_i(\mathbf{x}) = g(\mathbf{x})\} = 1
> $$

> 定理 1
> 
> 在满足假设1与2、且符合条件1与2的情况下，序列$\{\hat{\mathbf{x}}_k^\star\}_{k=0}^{\infty}$满足$\operatorname{Pr}\{\hat{\mathbf{x}}_{k}^{\star} \notin \mathcal{M}\text{ i.o.}\}=0$，其中$\mathcal{M} = \{\mathbf{x}: \mathbf{x}\in \Theta\text{ and either }\mathcal{N}(\mathbf{x})=\varnothing\text{ or }g(\mathbf{x}) \leq g(\mathbf{y})\text{ for all }\mathbf{y} \in \mathcal{N}(\mathbf{x})\}, \mathcal{N}(\mathbf{x}) = \{\mathbf{y}: \mathbf{y}\in \Theta\text{ and }\|\mathbf{x} - \mathbf{y}\| = 1\}$。

#### 证明思路

**第一步**，我们证明在符合条件1与2的情况下，当$\hat{\mathbf{x}}_{k}^{\star} = \mathbf{x}\text{ i.o.}$时，有$\operatorname{Pr}\{\mathcal{N}(\mathbf{x}) \not \subset \mathcal{E}(k)\text{ i.o. }\}=0$，即$\mathbf{x}$的所有“邻居”都会被仿真且估计。

这是由于$\operatorname{Pr}\{\mathbf{x}\in \mathcal{S}_k\}\geq \epsilon, \forall \mathbf{x} \in \mathcal{N}(\hat{\mathbf{x}}_k^\star) - \mathcal{E}(k - 1)$，因此$\lim\limits_{k\rightarrow \infty}\operatorname{Pr}\{\mathbf{y} \notin \mathcal{E}(k)\}=0, \forall \mathbf{y}\in \mathcal{N}(\hat{\mathbf{x}}_k^\star)$；又由于$\mathcal{N}(\hat{\mathbf{x}}_k^\star)$为有限集，故$\lim\limits_{k\rightarrow \infty}\operatorname{Pr}\{\mathcal{N}(\mathbf{x}) \not \subset \mathcal{E}(k)\}=0$；最后，由于$\operatorname{Pr}\{\mathcal{N}(\mathbf{x}) \not \subset \mathcal{E}(k)\text{ i.o. }\}\leq \operatorname{Pr}\{\mathcal{N}(\mathbf{x}) \not \subset \mathcal{E}(i)\}$，故第一步得证。

**第二步**，我们证明在满足假设1与2、且符合条件1与2的情况下，序列$\{\hat{\mathbf{x}}_k^\star\}_{k=0}^{\infty}$满足$\lim\limits_{k\rightarrow \infty} \left[g(\hat{\mathbf{x}}_{k}^{\star}) - \min_{\mathbf{y}\in \mathcal{E}_{k}} g(\mathbf{y})\right] = 0$ w.p. 1。

$$
\operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\min_{\mathbf{y}\in \mathcal{E}_{k}} g(\mathbf{y})\right| \geq \delta \text { i.o.}\right\}\\

\leq \operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\bar{G}_{k}\left(\hat{\mathbf{x}}_{k}^{\star}\right)\right| \geq \frac{\delta}{2} \text{ i.o.}\right\} + \operatorname{Pr}\left\{\left|\bar{G}_{k}\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\min_{\mathbf{y}\in \mathcal{E}_{k}} g(\mathbf{y})\right| \geq \frac{\delta}{2} \text{ i.o.}\right\}\\

\leq \operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\bar{G}_{k}\left(\hat{\mathbf{x}}_{k}^{\star}\right)\right| \geq \frac{\delta}{2} \text{ i.o.}\right\} + \operatorname{Pr}\left\{\left|\min _{\mathbf{y}\in \mathcal{E}_{k}} \bar{G}_{k}(\mathbf{y})-\min _{\mathbf{y}\in \mathcal{E}_{k}} g(\mathbf{y})\right| \geq \frac{\delta}{2} \text { i.o.}\right\}\\

\leq \operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\bar{G}_{k}\left(\hat{\mathbf{x}}_{k}^{\star}\right)\right| \geq \frac{\delta}{2} \text{ i.o.}\right\} + \operatorname{Pr}\left\{\max _{\mathbf{x} \in \mathcal{E}_{k}}\left|\bar{G}_{k}(\mathbf{x})-g(\mathbf{x})\right| \geq \frac{\delta}{2} \text { i.o.}\right\}\\

\leq 2 \operatorname{Pr}\left\{\mathcal{I}(\mathbf{x}\in \mathcal{E}_k)\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\delta}{2}\text{ for some }\mathbf{x} \in \mathcal{E}_{k} \text{ i.o. }\right\}\\

\leq 2 \operatorname{Pr}\left\{\mathcal{I}(\mathbf{x}\in \mathcal{E}_k)\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\delta}{2}\text{ for some }\mathbf{x} \in \mathcal{E}(\infty) \text{ i.o. }\right\}\\

\leq 2\sum\limits_{\mathbf{x}\in \mathcal{E}(\infty)}\operatorname{Pr}\left\{\mathcal{I}(\mathbf{x}\in \mathcal{E}_k)\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\delta}{2}\text{ i.o. }\right\}
$$

又由于条件1.1与2.2，再根据假设2，得到*（论文里的部分证明不是很理解为什么需要那么做）*
$$
\operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\min_{\mathbf{y}\in \mathcal{E}_{k}} g(\mathbf{y})\right| \geq \delta \text { i.o.}\right\}\\

\leq 2\sum\limits_{\mathbf{x}\in \mathcal{E}(\infty)}\operatorname{Pr}\left\{\mathcal{I}(\mathbf{x}\in \mathcal{E}_k)\left|g(\mathbf{x})-\bar{G}_{k}(\mathbf{x})\right| \geq \frac{\delta}{2}\text{ i.o. }\right\} = 0
$$

**第三步**，我们证明在满足假设1与2、且符合条件1与2的情况下，序列$\{g(\hat{\mathbf{x}}_k^\star)\}_{k=0}^{\infty}$收敛 w.p. 1。

根据第二步的结果、条件2.2和假设1，对于任意$\delta$，有
$$
\operatorname{Pr}\left\{g\left(\hat{\mathbf{x}}_{k}^{\star}\right) > g\left(\mathbf{x}_{0}\right) + \delta \text { i.o.}\right\}

\leq \operatorname{Pr}\left\{\left|g\left(\hat{\mathbf{x}}_{k}^{\star}\right)-\min_{\mathbf{y}\in \mathcal{E}_{k}} g(\mathbf{y})\right| \geq \delta \text { i.o.}\right\} = 0
$$

因此，
$$
\operatorname{Pr}\left\{\mathbf{x}_{k}^{\star} \notin \mathcal{L} \text { i.o.}\right\} = 0
$$

若$g(\mathbf{x}) = g(\mathbf{x}_0), \forall \mathbf{x}\in \mathcal{L}$，则$g\left(\hat{\mathbf{x}}_{k}^{\star}\right)$收敛至$g(\mathbf{x}_0)$w.p. 1；若$g(\mathbf{x}) \neq g(\mathbf{x}_0)$，则需要另外讨论。

同理，对于任意$\delta$，有
$$
\operatorname{Pr}\left\{g\left(\hat{\mathbf{x}}_{k}^{\star}\right) > g\left(\hat{\mathbf{x}}_{k-1}^{\star}\right) + \delta \text { i.o.}\right\} = 0
$$

令$\delta = \inf\{g(\mathbf{x})− g(\mathbf{y}): \mathbf{x},\mathbf{y} \in \mathcal{L}, g(\mathbf{x})− g(\mathbf{y}) > 0\}$，因此
$$
\operatorname{Pr}\left\{g\left(\hat{\mathbf{x}}_{k}^{\star}\right) > g\left(\hat{\mathbf{x}}_{k-1}^{\star}\right) \text { i.o.}\right\} = 0
$$

又根据假设1，$\mathcal{L}$是有限集，故序列$\{g(\hat{\mathbf{x}}_k^\star)\}_{k=0}^{\infty}$收敛 w.p. 1。

**第四步**，证明定理1。根据第三步，我们得到
$$
\operatorname{Pr}\left\{\mathbf{x}_{k}^{\star} \notin \mathcal{L} \text { i.o.}\right\} = 0
$$

又根据假设1，$\mathcal{L}$是有限集，因此一定能找到$\mathbf{x}$，使得$\hat{\mathbf{x}}_{k}^{\star} = \mathbf{x}$ i.o. ，记$\hat{\mathbf{x}}_{k_i}^{\star} = \mathbf{x}, i = 1, 2, \cdots$，于是有$\lim\limits_{i\rightarrow \infty} \left[g(\hat{\mathbf{x}}_{k_i}^{\star}) - \min_{\mathbf{y}\in \mathcal{E}_{k_i}} g(\mathbf{y})\right] = 0$ w.p. 1且$\lim\limits_{i\rightarrow \infty} g(\hat{\mathbf{x}}_{k}^{\star}) = g(\mathbf{x})$，所以
$$
\lim\limits_{i\rightarrow \infty} \left[\min_{\mathbf{y}\in \mathcal{E}_{k_i}} g(\mathbf{y})\right] = g(\mathbf{x})\text{ w.p. 1 }
$$

又根据第一步的结论，$\operatorname{Pr}\{\mathcal{N}(\mathbf{x}) \not \subset \mathcal{E}(k)\text{ i.o. }\}=0$，又由于条件2.2，有$\operatorname{Pr}\{\mathcal{N}(\mathbf{x}) \not \subset \mathcal{E}_{k_i}\text{ i.o. }\}$，因此
$$
\min_{\mathbf{y}\in \mathcal{N}(\mathbf{x})} g(\mathbf{y})\geq \lim\limits_{i\rightarrow \infty} \left[\min_{\mathbf{y}\in \mathcal{E}_{k_i}} g(\mathbf{y})\right] = g(\mathbf{x})\text{ w.p. 1 }
$$

### 改进COMPASS算法

根据LCRS算法的框架，在估计阶段（estimation scheme），无需对$\mathcal{V}_k$中所有的解进行估计，只需要对$\mathbf{x}_{0}, \mathcal{N}\left(\widehat{\mathbf{x}}_{k-1}^{*}\right) \cap \mathcal{E}(k-1), \mathcal{S}_{k}$这三个集合中的解进行估计即可。*（这里不是很理解文章中为何要构造$\mathcal{P^\prime}$，计算机中对$\mathcal{N}\left(\widehat{\mathbf{x}}_{k-1}^{*}\right) \cap \mathcal{E}(k-1)$的查找可能不难实现？也有可能是使算法更加高效？）*

## AHA算法

对于完全约束问题，LCRS的条件可以进行调整，由此可以产生更高效的AHA算法。

#### 条件

算法需要满足的条件调整为如下，其中，条件1比原LCRS更强，而条件2更弱：

> 条件 1
> 
> - 抽样分布$F_k$满足$\operatorname{Pr}\{\mathbf{x}\in \mathcal{S}_k\}\geq \epsilon, \forall \mathbf{x} \in \mathcal{N}(\hat{\mathbf{x}}_k^\star)$，其中$\epsilon$与$k$无关。

> 条件 2
>
> - $\mathcal{E}_{k}\subset \mathcal{S}(k)$；
> - $\mathcal{E}_{k}$包括$\mathbf{x}_{k-1}^{*}, \mathcal{S}_{k}$；
> - $a_{k}(\mathbf{x})$满足$\min _{\mathbf{x} \in \mathcal{E}_{k}} N_{k}(\mathbf{x}) \geq 1, \forall k=1,2, \cdots$且$\min _{\mathbf{x} \in \mathcal{E}_{k}} N_{k}(\mathbf{x}) \rightarrow \infty$ w.p. 1 as $k \rightarrow \infty$，其中$N_{k}(\mathbf{x})$为$k$次迭代共计在解$\mathbf{x}$处进行仿真的次数。

AHA算法的MPA和COMPASS的MPA有所不同，AHA的MPA是解空间各个维度上距离当前找到的最好的解最近的已探索过的解包围的超盒，即：
> For a visited solution $\mathbf{x}$, let $x^{(d)}$ be its $d$ th coordinate, $1 \leq d \leq D$. Let $l_{k}^{(d)}=\max _{\mathrm{x} \in \mathcal{S}(k), \mathrm{x} \neq{\hat{x}}_{k}^{*}}\left\{x^{(d)}: x^{(d)}<\hat{x}^{*(d)}\right\}$ if it exists; otherwise, let $l_{k}^{(d)}=-\infty$. Similarly, let $u_{k}^{(d)}=$ $\min _{\mathbf{x} \in \mathcal{S}(k), \mathbf{x} \neq \hat{x}_{k}^{*}}\left\{x^{(d)}: x^{(d)}>\hat{x}_{k}^{*(d)}\right\}$ if it exists; otherwise, let $u_{k}^{(d)}=\infty$. Note that $u_{k}^{(d)}$ and $l_{k}^{(d)}$ may be $\pm \infty$ because $\hat{\mathbf{x}}^{*}$ may be on the boundary, or we may not have visited enough solutions yet. The hyperbox containing $\hat{\mathbf{x}}_{k}^{*}$ is $\mathcal{H}_{k}=\left\{\mathbf{x}: l_{k}^{(d)} \leq x^{(d)} \leq u_{k}^{(d)}, 1 \leq d \leq D\right\}$.

最初的COMPASS算法和AHA算法构造的最有前景区域之间的比较如下图所示。

<img src="docs/COMPASS/aha.png" alt="最有前景区域" width="400"/>

## 仿真优化解决港口泊位分配问题

### 问题背景

停泊在港口的船只包括大船（deep-sea vessel）与小船（feeder）。大船指在深海运输集装箱的船只，大船的到达时间点与装（卸）货需要花费的时间是确定的，大船会占据几个泊位；小船指在近海运输货物的船只，每只小船只有在一定时间段内才被允许到达，它的装（卸）货时间在小船到达前是未知的（服从一定的概率分布）、在小船到达后就是确定的，一艘小船占据一个泊位。大船到其停泊时间后，则进入对应泊位停泊；小船遵循先到先进入的原则，选择空泊位停泊（若直到小船离开前的时间里所在泊位已被大船“预订”，则该泊位也不能停泊）。

- 决策内容包括：
  - 大船停泊位置与停泊时间；
  - 小船到达时间。
- 约束内容包括：
  - 必须为大船安排泊位与停泊时间；
  - 必须为小船安排到达时间；
  - 大船的泊位与停泊安排不能冲突；
  - 期望排队的小船数不能超过一定值。
- 大船有计划离开时间，小船有计划到达时间。最终的目标为大船的实际离开时间不能过晚、小船的实际到达时间不能偏差计划过大。两者加权，得到目标函数。

### 问题建模

$$
\mathbf{P}: \operatorname{minimize} \sum_{i=1}^{N_{1}} C_{1 i} l_{1 i}+\sum_{i=1}^{N_{2}} C_{2 i} l_{2 i} \\
\text { subject to } \sum_{b=1}^{B-R_{i}+1} \sum_{t=A_{i}}^{T-H_{i}} x_{i b t}=1 \quad\left(i=1, \ldots, N_{1}\right) \\
\sum_{t=\underline{E_{i}}}^{\overline{{E}_{i}}} y_{i t}=1 \quad\left(i=1, \ldots, N_{2}\right) \\
\sum_{i=1}^{N_{1}} \sum_{b^{\prime}=\max \left\{1, b-R_{i}+1\right\}}^{\min \left\{b, B-R_i+1\right\}} \sum_{t^{\prime}=\max \left\{0, t-H_{i}+1\right\}}^{\min \left\{t, T-H_i\right\}} x_{i b^{\prime} t^{\prime}} \leq 1 \quad \text{for any }(b, t) \\
E\left[Q_{t}(\mathbf{x}, \mathbf{y})\right] \leq \overline{Q} \quad(t=0,1, \ldots, T-1)\\
l_{1 i}=\max \left\{0, \sum_{b=1}^{B-R_{i}+1} \sum_{t=A_{i}}^{T-H_{i}}\left(t+H_{i}\right) x_{i b t}-D_{i}\right\} \quad\left(i=1, \ldots, N_{1}\right) \\
    l_{2 i}=\left|\sum_{t=\underline{E_{i}}}^{\overline{{E}_{i}}} t y_{i t}-S_{i}\right| \quad\left(i=1, \ldots, N_{2}\right) \\
    x_{i b t} \in\{0,1\} \quad \text{for any }(i, b, t) \\
    y_{i t} \in\{0,1\} \quad \text{for any }(i, t)
$$

这是一个**随机规划**模型，其中
- 即使$N_2 = 0$，该问题也是一个NP-难问题。
- $E\left[Q_{t}(\mathbf{x}, \mathbf{y})\right]$不能显式表达（除非极特殊的情况）；

因此需要设计仿真优化算法进行求解。

### 问题求解

#### 对原问题进行松弛

化原来的“硬”约束为“软”约束，去除约束（5），将目标变为$\text { minimize } \sum_{i=1}^{N_{1}} C_{1 i} l_{1 i}+\sum_{i=1}^{N_{2}} C_{2 i} l_{2 i}+\lambda \Delta$，其中$\Delta=\max \left\{0, \max _{0,1, \ldots, T-1}\left\{E\left[Q_{t}(\mathbf{x}, \mathbf{y})\right]\right\}-\overline{Q}\right\}$，新的优化问题记作$\mathbf{P}^{\prime}(\lambda)$。

#### 三阶段仿真优化方法

1. 全局阶段：利用自适应遗传算法得到松弛问题的若干较优解；
2. 局部阶段：将较优解聚成$K$簇，利用自适应超盒算法（Adaptive hyperbox algorithm）寻找每簇内的较优解，其中“超盒”即聚成的“簇”，也就是全局阶段寻找到的较优解的邻域；
3. 清理阶段：对于局部阶段求得的解应用排序-选择方法（Ranking and Selection）得到最优解。

#### 全局阶段——自适应遗传算法

- 编码
  - $N_1 + N_2$的全排列$\mathbf{p}=\left(p_{1}, \ldots, p_{N_{1}+N_{2}}\right)$。其中$N_1$为大船数，$N_2$为小船数。
- 解码
  - 如果$p_{i} \leq N_{1}$，则代表着大船$p_{i}$；否则，则代表着小船$p_{i} - N_{2}$；
  - 从$p_{1}$到$p_{N_1 + N_2}$，按此优先级顺序，贪婪地安排大船的泊位、停泊时间和小船的到达时间；
  - 安排方法：满足一定要求的情况下（大船的到达时间、小船的到达时间段），有空就填、尽前（即时间尽可能靠前、泊位号尽可能小）。
- 初始化
  - 设置$\lambda = \lambda_0$，生成初始种群，包含$s_1$个个体。生成每一个个体时，对于每一个$i=1, \ldots, N_{1}+N_{2}$，随机抽取染色体中的空位置放入。
- 交叉
  - 选择$s_2$对染色体进行交叉，对于每一对染色体$\mathbf{p}=\left(p_{1}, \ldots, p_{N_{1}+N_{2}}\right)$和$\overline{\mathbf{p}}=\left(\overline{p}_{1}, \ldots, \overline{p}_{N_{1}+N_{2}}\right)$，随机生成两个位点$j$和$j^{\prime}$，记$\mathbf{q} = (\underbrace{0, \ldots, 0}_{j-1 \text { times }}, p_{j}, p_{j+1}, \ldots, p_{j^{\prime}}, \underbrace{0, \ldots, 0}_{N_{1}+N_{2}-j^{\prime} \text { times }}),  \overline{\mathbf{q}} = (\underbrace{0, \ldots, 0}_{j-1 \text { times }}, \overline{p}_{j}, \overline{p}_{j+1}, \ldots, \overline{p}_{j^{\prime}}, \underbrace{0, \ldots, 0}_{N_{1}+N_{2}-j^{\prime} {\text { times }}})$，将$\mathbf{q}$中没有的数字按$\overline{\mathbf{p}}$中这些数字的顺序填入，将$\overline{\mathbf{q}}$中没有的数字按$\mathbf{p}$中这些数字的顺序填入。
- 变异
  - 对于群体中的每一条染色体，以一定概率交换随机的两个位点。
- 计算适应度
  - 对于群体中每一条染色体，解码成问题$\mathbf{P}^{\prime}(\lambda)$的解，进行一定次数（$n_1$）的仿真，适应度即为众多仿真结果的平均值。
- 选择
  - 保留适应度高的$s_1$条染色体。
- 更新$\lambda$
  - 如果对于适应度最高的染色体来说，有$\max _{0,1, \ldots, T-1}\left\{E\left[Q_{t}(\mathbf{x}(\mathbf{p}), \mathbf{y}(\mathbf{p}))\right]\right\}$，则更新$\lambda = (1+\zeta) \lambda$。

> 定理 1
> 
> 如果原问题$\mathbf{P}$是可行的，则一定存在$\hat{\lambda}>0$，使得问题$\mathbf{P}^{\prime}(\hat{\lambda})$的最优解就是原问题的最优解。

#### 局部阶段——自适应超盒算法

令$e_{1 i}(\mathbf{x})$表示大船$i$占据的第一个泊位序号，$e_{2 i}(\mathbf{x})$表示为大船$i$安排的停泊时间点，$e_{3 i}(\mathbf{y})$表示为小船$i$安排的到达时间点。

定义两个解之间的距离$d_{j j^{\prime}}$为
$$
\sqrt{\sum_{i=1}^{N_{1}}\left[\frac{e_{1 i}\left(\mathbf{x}_{j}\right)-e_{1 i}\left(\mathbf{x}_{j^{\prime}}\right)}{B}\right]^{2}+\sum_{i=1}^{N_{1}}\left[\frac{e_{2 i}\left(\mathbf{x}_{j}\right)-e_{2 i}\left(\mathbf{x}_{j^{\prime}}\right)}{T}\right]^{2}+\sum_{i=1}^{N_{2}}\left[\frac{e_{3 i}\left(\mathbf{y}_{j}\right)-e_{3 i}\left(\mathbf{y}_{j^{\prime}}\right)}{T}\right]^2}
$$

以$(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)$作为变量进行AHA算法实现，一方面，数据降维，变量的维度从$(N_1B + N_2)T$降至$2N_1 + N_2$，且目标函数更加“平滑”；另一方面，在超盒中随机抽样到的解可能不满足约束：$\sum_{i=1}^{N_{1}} \sum_{b^{\prime}=\max \left\{1, b-R_{i}+1\right\}}^{\min \left\{b, B-R_i+1\right\}} \sum_{t^{\prime}=\max \left\{0, t-H_{i}+1\right\}}^{\min \left\{t, T-H_i\right\}} x_{i b^{\prime} t^{\prime}} \leq 1$。若在超盒中随机抽样到的解可能不满足该约束，则对该解进行“微调”，改变安排冲突的大船的停泊时间为$\bar{t}_{i} \leftarrow \min \left\{t \mid \sum_{b^{\prime}=\bar{b}_{i}}^{\bar{b}_{i}+R_{i}-1} \sum_{t^{\prime}=t}^{t+H_{i}-1} \omega_{b^{\prime} t^{\prime}}=0 ; t=A_{i}, A_{i}+1, \ldots, T-H_{i}+1\right\}$，其中$\omega_{bt}$为$t$时间占据泊位$b$的大船数。但*如此进行微调是否还能保证“平滑”*？

## 总结

- Jia等首先通过引入惩罚系数，化“硬”约束为“软”约束对原问题进行了**松弛**，同时随着迭代过程的进行，惩罚系数的值会动态递增，有利于问题的求解的同时可以使得最后的解尽可能满足约束。
- Jia等有效地实现了变量的**降维**。在全局阶段，用$N_1 + N_2$的全排列作为遗传算法的编码；在局部阶段，以$(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)$这个$2N_1 + N_2$维向量作为变量进行聚类和自适应超盒算法实现。从而，一方面有利于加速迭代过程，另一方面使得目标函数更加“平滑”。
- Hong等对最初的COMPASS算法不断进行改进，主要体现在：
  - Hong等通过LCRS算法框架对COMPASS算法的估计阶段进行了改进:一方面，不需要对**所有**抽样抽到的解进行大量仿真；另一方面，寻找到近似有效约束点集，只需要对该集合中的点、新抽样的点和初始点进行多次仿真即可保证局部收敛性；
  - Xu等通过求解线性规划问题得到了精确的有效约束点集，仿真耗费的时间大于求解小型线性规划，因此该步骤可以有效提升效率；
  - Xu等在解集有限的情况下，调整了LCRS的条件，并调整了最有前景区域，产生了更高效的自适应超盒算法。

- Jia等在局部阶段，以$(\mathbf{e}_1, \mathbf{e}_2, \mathbf{e}_3)$作为变量进行自适应超盒算法实现，会导致大船的停泊安排发生冲突，作者通过把冲突大船的停泊时间安排到可行的最前的时间进行“微调”。但是，我们认为这种“微调”会导致目标函数不再“平滑”，可以考虑把冲突大船的停泊时间安排到距离当前冲突时间最近的可行时间。
- Jia等在局部阶段使用的自适应超盒算法仍然是对过去迭代**所有**抽样抽到的解进行仿真，原因不明。我们认为可以考虑依据Xu等结论，仅对上一轮迭代后寻找到的当前最优解和新抽样的解进行仿真。
- Hong等依据LCRS算法框架对COMPASS算法进行了改进，但改进COMPASS算法中对估计阶段所要进行仿真的解是多于LCRS算法框架中局部收敛条件要求的解的，我们不知道这种做法的目的。我们认为可以考虑把仅对LCRS算法框架中局部收敛条件要求的解进行仿真和目前的改进COMPASS算法进行实验对比。

## 参考文献

- Boesel, J., Nelson, B.L., Kim, S.-H., 2003. Using Ranking and Selection to “Clean Up” after Simulation Optimization. Operations Research 51, 814–825.
  - https://doi.org/10.1287/opre.51.5.814.16751
- **Jeff Hong, L., Nelson, B.L., 2006. Discrete optimization via simulation using COMPASS. Operations Research 54, 115–129.**
  - https://doi.org/10.1287/opre.1050.0237
- **Jeff Hong, L., Nelson, B.L., 2007. A framework for locally convergent random-search algorithms for discrete optimization via simulation. ACM Trans. Model. Comput. Simul. 17, 19.**
  - https://doi.org/10.1145/1276927.1276932
- **Jia, S., Li, C.-L., Xu, Z., 2020. A simulation optimization method for deep-sea vessel berth planning and feeder arrival scheduling at a container port. Transportation Research Part B: Methodological 142, 174–196.**
  - https://doi.org/10.1016/j.trb.2020.10.007
- Pichitlamken, J., Nelson, B.L., 2002. A combined procedure for optimization via simulation, in: Proceedings of the Winter Simulation Conference. Presented at the Proceedings of the Winter Simulation Conference, pp. 292–300 vol.1.
  - https://doi.org/10.1109/WSC.2002.1172898
- Xu, J., Nelson, B.L., Hong, J.L., 2010. Industrial strength COMPASS: A comprehensive algorithm and software for optimization via simulation. ACM Trans. Model. Comput. Simul. 20, 1–29.
  - https://doi.org/10.1145/1667072.1667075
- **Xu, J., Nelson, B.L., Hong, L.J., 2013. An Adaptive Hyperbox Algorithm for High-Dimensional Discrete Optimization via Simulation Problems. INFORMS Journal on Computing 25, 133–146.**
  - https://doi.org/10.1287/ijoc.1110.0481
