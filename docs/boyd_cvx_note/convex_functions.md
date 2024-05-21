# Convex Functions

## Basic properties and examples

### Definition

A function $f:\mathbb{R}^n\rightarrow\mathbb{R}$ is **convex** if $\operatorname{dom}f$ is a *convex set* and if for all $x, y \in \operatorname{dom}f$, and $\theta$ with $0 \leq \theta \leq 1$, we have
$$
f(\theta x + (1 - \theta)y) \leq \theta f(x) + (1 - \theta)f(y)
$$
A function $f$ is **strictly convex** if strict inequality holds whenever $x \neq y$ and $0 < \theta < 1$. $f$ is **concave** if $-f$ is convex, and strictly concave if $-f$ is strictly convex.

- Methods to check convexity:

  1. Definition;
  2. $f$ is convex if and only if $g(t) = f(x+tv), x \in \operatorname{dom}f, x+tv\in \operatorname{dom}f$;
  3. First-order conditions: Suppose $f$ is differentiable, then $f$ is convex if and only if $\operatorname{dom} f$ is convex and $\forall x,y \in \operatorname{dom}f, f(y) \geq f(x) + \nabla f(x)^T (y - x)$;
  4. Second-order conditions: Suppose $f$ is twice differentiable, then $f$ is convex if and only if $\operatorname{dom} f$ is convex and $\forall x \in \operatorname{dom}f, \nabla^2 f(x) \geq 0$;
     - Example: *Quadratic functions* $f(x) = (1/2)x^T Px + q^T x + r$;
  5. Epigraph: A function is convex if and only if its epigraph $\operatorname{epi}f=\{(x,t)|x\in \operatorname{dom}f, f(x)\leq t\}$, is a convex set;
     - Example: $f(x, Y) = x^T Y^{-1} x, \operatorname{dom}f = \mathbb{R}^n\times\mathbb{S}^n_{++}$;
     - The relationship between the first order condition $f(y) \geq f(x) + \nabla f(x)^T (y - x)$ and the supporting hyperplane at $(x, f(x))$($(y, t) \in \operatorname{epi} f \Rightarrow\left[\begin{array}{c}\nabla f(x) \\ -1\end{array}\right]^{T}\left(\left[\begin{array}{c}y \\ t\end{array}\right]-\left[\begin{array}{c}x \\ f(x)\end{array}\right]\right) \leq 0$).
     - Sublevel sets: $C_{\alpha} =\{x\in \operatorname{dom} f \mid f(x)\leq\alpha\}$. Sublevel sets of a convex function are convex (often use the property to check that a function is *not* convex).

### Examples

Examples below are all convex functions:
- Exponential: $f(x) = e^{ax}, x\in \mathbb{R}$;
- Powers: $f(x) = x^{a}, x\in \mathbb{R}_{++}, a \geq 1 \text{ or } a \leq 0$;
- Powers of absolute value: $f(x) = |x|^{a}, x\in \mathbb{R}, a \geq 1$;
- Logarithm: $f(x) = -\log x, x\in \mathbb{R}_{++}$;
- Negative entropy$f(x) = -x\log x, x\in \mathbb{R}_{+}$;
- Norms(proved by definition);
- Max function: $f(x) = \max_i x_i$(proved by definition);
- Quadratic-over-linear function: $f(x, Y) = x^T Y^{-1} x, \operatorname{dom}f = \mathbb{R}^n\times\mathbb{S}^n_{++}$(proved by second-order condition);
- Log-sum-exp: $f(x) = \log (e^{x_1} + \cdots + e^{x_n})$(proved by second-order condition and Cauchy-Schwarz inequality);
- Negative geometric mean: $f(x) = -(\prod_{i = 1}^n x_i)^{1/n}, x\in \mathbb{R}_{++}^n$(proved by second-order condition and Cauchy-Schwarz inequality);
- Log-determinant: $f(x) = -\log \det X, X\in \mathbb{S}_{++}^n$;
  - $$
    g(t) = f(X + tV) = -\log \det (X + tV)\\
    = -\log \det X^{1/2}(I + tX^{-1/2}VX^{-1/2})X^{1/2} = -\log \det X - \log \det (I + tX^{-1/2}VX^{-1/2})\\
    = -\log \det X -\log \det Q(I + t\Lambda)Q^{-1} = -\log \det X -\sum_{i=1}^n\log \det(1 + t\lambda_i)
    $$
  - Apply second-order condition on $g(t)$.

### Jensen’s inequality

If $x$ is a random variable such that $x \in \operatorname{dom}f$ with probability one, and $f$ is convex, then we have
$$
f(\operatorname{E}x) ≤ \operatorname{E} f(x)
$$
provided the expectations exist. Many famous inequalities can be derived by applying Jensen’s inequality to some appropriate convex function.

## Operations that preserve convexity

Operations that preserve convexity include:
- Nonnegative weighted sums: $f = w_1f_1 + \cdots + w_mf_m$, $f_i$ is convex;
- Composition with an affine mapping: $g(x) = f(Ax + b)$, $f$ is convex;
- Pointwise maximum: $f(x) = \max\{f_1(x),\cdots,f_m(x)\}$, $f_i$ is convex;
  - Piecewise-linear functions $f(x)=\max \left\{a_{1}^{T} x+b_{1}, \cdots, a_{L}^{T} x+b_{L}\right\}$;
  - Sum of $r$ largest components of a vector;
- Pointwise supremum: $g(x) = \sup\limits_{y\in \mathcal{A}} f(x, y)$, $f(x,y)$ is convex in $x$ (can be proved by epigraph);
  - Support function of a set $S_C(x) = \sup\{x^T y | y \in C\}$;
  - Distance to farthest point of a set $f(x) = \sup\limits_{y\in \mathcal{A}} \Vert x - y \Vert$;
  - Maximum eigenvalue of a symmetric matrix $f(X)=\sup \left\{y^{T} X y \mid\|y\|_{2}=1\right\}$;
- Minimization: If $f$ is convex in $(x,y)$, and $C$ is a convex nonempty set, then the function $g(x) = \inf\limits_{y\in C} f(x,y)$;
  - Schur complement: Suppose the quadratic function $f(x, y)=x^{T} A x+2 x^{T} B y+y^{T} C y$ (where $A$ and $C$ are symmetric) is convex in $(x, y)$, which means
    $$
    \left[\begin{array}{cc}
    A & B \\
    B^{T} & C
    \end{array}\right] \succeq 0 .
    $$
    We can express $g(x)=\inf _{y} f(x, y)$ as
    $$
    g(x)=x^{T}\left(A-B C^{\dagger} B^{T}\right) x
    $$
    By the minimization rule, $g$ is convex, so we conclude that $A-B C^{\dagger} B^{T} \succeq 0$. If $C$ is invertible, i.e., $C \succ 0$, then the matrix $A-B C^{-1} B^{T}$ is called the Schur complement of $C$ in the matrix
    $$
    \left[\begin{array}{cc}
    A & B \\
    B^{T} & C
    \end{array}\right]
    $$
  - Distance to a set: $f(x) = \inf\limits_{y\in S} \Vert x-y\Vert$;
- Composition under some conditions;
- Perspective of a function: If $f$ is a convex function, then so is its perspective function $g(x, t) = tf(x/t)$.