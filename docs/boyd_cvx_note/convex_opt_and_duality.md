# Convex Optimization and Duality

## Convex Optimization

### Standard form of convex optimization

A convex optimization problem is one of the form:
$$
\begin{array}{ll}
\operatorname{minimize} & f_{0}(x) \\
\text {subject to } & f_{i}(x) \leq 0,  i=1, \cdots, m \\
& a_{i}^{T} x=b_{i},  i=1, \cdots, p
\end{array}
$$
where $f_0, \cdots, f_m$ are convex functions. The **feasible set** of a convex optimization problem is convex. The requirements of a convex optimization problem are:

- the objective function must be convex;
- the inequality constraint functions must be convex;
- the equality constraint functions must be affine.

### Equivalent convex problems

- Eliminating equality constraints;
  - $Ax = b \Leftrightarrow x = Fz + x_0$;
- Introducing equality constraints;
- Slack variables;
- Epigraph problem form;
  - $$
    \begin{array}{ll}\operatorname{minimize} & t \\
    \text {subject to } & f_{0}(x)-t \leq 0\end{array}
    $$
- Minimizing over some variables.

### Local and global optima

Any locally optimal point of a convex optimization problem is also (globally) optimal. (proof by contradiction)

### Optimality criterion

Suppose that the objective $f_0$ in a convex optimization problem is differentiable, then $x$ is optimal if and only if $x$ is feasible and
$$
\nabla f_{0}(x)^{T}(y-x) \geq 0
$$
for all $y \in X$.

Simple examples of optimality criterion:
- Unconstrained problems: $\nabla f_{0}(x) = 0$;
- Problems with equality constraints $Ax - b = 0$ only: $y = x + v, v \in \mathcal{N}(A)$, then $\nabla f_{0}(x)^{T}v \geq 0, \forall v \in \mathcal{N}(A)$, so we get $\nabla f_{0}(x)^{T}v = 0, \forall v \in \mathcal{N}(A) \Leftrightarrow \nabla f_{0}(x)\in \mathcal{R}(A^T)\Leftrightarrow \exist v\in \mathbb{R}^p, \nabla f_{0}(x) + A^Tv = 0$;
- Minimization over the nonnegative orthant (i.e. $x\geq 0$): $\nabla f_{0}(x)^{T}(y-x) \geq 0, \forall y \geq 0$, so $\nabla f_{0}(x)\geq 0, -\nabla f_{0}(x)^{T}x \geq 0$. But $x\geq 0$, we finally get that $\nabla f_{0}(x)\geq 0, (\nabla f_{0}(x))_ix_i = 0$.
- We will discuss more generic optimality conditions later.

## The Lagrange dual function

We consider an optimization problem in the standard form:
$$
\begin{array}{ll}
\operatorname{minimize} & f_{0}(x) \\
\text {subject to } & f_{i}(x) \leq 0,  i=1, \cdots, m \\
& h_{i}(x)=0,  i=1, \cdots, p
\end{array}
$$

**Lagrangian**: $L(x, \lambda, \nu)=f_{0}(x)+\sum\limits_{i=1}^{m} \lambda_{i} f_{i}(x)+\sum\limits_{i=1}^{p} \nu_{i} h_{i}(x)$.

**Lagrange dual function**: $g(\lambda, \nu)=\inf\limits_{x \in \mathcal{D}} L(x, \lambda, \nu)=\inf\limits_{x \in \mathcal{D}}\left(f_{0}(x)+\sum\limits_{i=1}^{m} \lambda_{i} f_{i}(x)+\sum\limits_{i=1}^{p} \nu_{i} h_{i}(x)\right)$.
  - Concave (the pointwise infimum of a family of affine functions of $(\lambda,\nu)$);
  - For any $\lambda \geq 0$ and any $\nu$ we have $g(\lambda, \nu) \leq p^{\star}$.

### Examples:
- $\min x^Tx, \text{ s.t. } Ax = b$;
  - $L(x, \nu) = x^Tx + \nu^T(Ax - b)$;
  - $g(\nu) = L(-(1/2)A^T\nu, \nu) = -(1/4)\nu^TAA^T\nu - b^T\nu$;
- $\min c^Tx, \text{ s.t. } Ax = b, x\geq 0$;
  - $L(x, \lambda, \nu) = c^Tx - \lambda^T x + \nu^T(Ax - b) = (c- \lambda + A^T\nu)^Tx - b^T\nu$;
  - $g(\lambda, \nu) = \begin{cases}
      -b^T\nu, A^T\nu - \lambda + c = 0\\
      -\infty, \text{ otherwise }
    \end{cases}$;
- $\min x^TWx, \text{ s.t. } x_i^2 = 1$;
  - $L(x, \nu) = x^TWx + x^T\operatorname{diag}(\nu)x - 1^T\nu = x^T(W + \operatorname{diag}(\nu))x - 1^T\nu$;;
  - $g(\nu) = \begin{cases}
      -1^T\nu, W + \operatorname{diag}(\nu) \geq 0\\
      -\infty, \text{ otherwise }
    \end{cases}$;

## The Lagrange dual problem

Since for any $\lambda \geq 0$ and any $\nu$ we have $g(\lambda, \nu) \leq p^{\star}$, we can get a good lower bound of the **primal problem** by solving:
$$
\begin{array}{ll}
\text{maximize } & g(\lambda, \nu)\\
\text{subject to } & λ \geq 0
\end{array}
$$
which is called the **Lagrange dual problem** associated with the primal problem.

- **Dual feasible**: $(\lambda, \nu)$ is feasible for the dual problem, that is  $\lambda \geq 0$ and $g(\lambda, \nu) > −\infty$;

- Example of LP (linear programming): $\min c^Tx, \text{ s.t. } Ax = b, x\geq 0$;
  - $L(x, \lambda, \nu) = c^Tx - \lambda^T x + \nu^T(Ax - b) = (c- \lambda + A^T\nu)^Tx - b^T\nu$;
  - $g(\lambda, \nu) = \begin{cases}
      -b^T\nu, A^T\nu - \lambda + c = 0\\
      -\infty, \text{ otherwise }
    \end{cases}$;
  - Dual problem: 
    $$
    \begin{array}{ll}
    \text{maximize } & g(\lambda, \nu)\\
    \text{subject to } & λ \geq 0
    \end{array}\\
    \Leftrightarrow
    \begin{array}{ll}
    \text{maximize } & -b^T\nu\\
    \text{subject to } & \begin{cases}
        λ \geq 0\\
        A^T\nu - \lambda + c = 0
    \end{cases}\Leftrightarrow A^T\nu + c \geq 0
    \end{array}
    $$
  - Weak duality: $d^{\star} \leq p^{\star}$, where $d^{\star}$ is the optimal value of the Lagrange dual problem and $p^{\star}$ the primal problem.

### Strong duality and Slater’s constraint qualification

**Strong duality**: $d^{\star} = p^{\star}$. Strong duality does *not*, in general, hold.

**Slater’s condition**: There exists an $x \in \operatorname{relint} \mathcal{D}$ such that
$$
f_{i}(x)<0,  i=1, \cdots, m,  A x=b
$$
- Such a point is sometimes called **strictly** feasible;
- Slater's theorem states that strong duality holds if:
  1. Slater's condition holds;
  2. the problem is convex;
- If the first $k$ constraint functions $f_{1}, \cdots, f_{k}$ are affine, then strong duality holds provided the following weaker condition holds: There exists an $x \in \operatorname{relint} \mathcal{D}$ such that $f_{i}(x) \leq 0, i=1, \cdots, k; f_{i}(x)<0, i=k+1, \cdots, m; Ax=b$.

## Optimality conditions

For ***any*** optimization problem with differentiable objective and
constraint functions for which ***strong duality obtains***, any pair of primal and dual optimal points must satisfy the **KKT conditions**:
$$
f_{i}\left(x^{\star}\right) \leq 0, i=1, \cdots, m \\
h_{i}\left(x^{\star}\right) = 0, i=1, \cdots, p \\
\lambda_{i}^{\star}  \geq 0, i=1, \cdots, m \\
\lambda_{i}^{\star} f_{i}\left(x^{\star}\right) = 0, i=1, \cdots, m \\
\nabla f_{0}\left(x^{\star}\right)+\sum_{i=1}^{m} \lambda_{i}^{\star} \nabla f_{i}\left(x^{\star}\right)+\sum_{i=1}^{p} \nu_{i}^{\star} \nabla h_{i}\left(x^{\star}\right) = 0
$$

For ***any convex*** optimization problem with differentiable objective and constraint functions, any points that satisfy the KKT conditions are primal and dual optimal.

- Complementary slackness (under strong duality condition): 
  - $$
    \begin{aligned}
    f_{0}\left(x^{\star}\right) &=g\left(\lambda^{\star}, \nu^{\star}\right) \\
    &=\inf _{x}\left(f_{0}(x)+\sum_{i=1}^{m} \lambda_{i}^{\star} f_{i}(x)+\sum_{i=1}^{p} \nu_{i}^{\star} h_{i}(x)\right) \\
    & \leq f_{0}\left(x^{\star}\right)+\sum_{i=1}^{m} \lambda_{i}^{\star} f_{i}\left(x^{\star}\right)+\sum_{i=1}^{p} \nu_{i}^{\star} h_{i}\left(x^{\star}\right) \\
    & \leq f_{0}\left(x^{\star}\right)
    \end{aligned}\\
    \Rightarrow \sum_{i=1}^{m} \lambda_{i}^{\star} f_{i}\left(x^{\star}\right)=0
    $$