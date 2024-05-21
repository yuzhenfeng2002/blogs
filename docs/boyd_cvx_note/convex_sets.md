# Convex Sets

## Affine and convex sets

**Affine Set**: $C \subseteq \mathbb{R}^n$ is an affine set if and only if $\forall x_1, x_2 \in C, \theta \in \mathbb{R}, \theta x_1 + (1 - \theta) x_2 \in C$.

- Solution set of a linear equation(i.e. $C = \{x|Ax = b\}$) is an affine set.

**Convex Set**: $C \subseteq \mathbb{R}^n$ is an affine set if and only if $\forall x_1, x_2 \in C, \theta \in [0, 1], \theta x_1 + (1 - \theta) x_2 \in C$.

- An affine set is a convex set.

We call a point of the form $\theta_1x_1 + \cdots + \theta_kx_k$, where $\theta_1 + \cdots + \theta_k = 1, \theta_i \geq 0$ a *convex combination* of points $x_1, \cdots, x_k$. The *convex hull* of a set C, denoted $\mathbf{conv} C$ is:
$$
\mathbf{conv} C = \{\theta_1x_1 + \cdots + \theta_kx_k|x_i \in C, \theta_i \geq 0, \theta_1 + \cdots + \theta_k = 1\}
$$

## Some important examples

- Affine sets(hence convex)
  
  - the empty set, a single point, line, subspace, the whole space
  - hyperplane: $\{x|a^T x = b\}$

- Other convex sets
  
  - line segment, a ray
  - halfspace: $\{x|a^T x \leq b\}$
  - ellopsoid: $\{x|(x - x_c)^TP^{-1}(x - x_c) \leq 1\}, P \in \mathbb{S}^n_{++}$
    - when $P = E$, it is an euclidean ball
    - another common representation of an ellipsoid is $\{x_c + Au| \left\Vert u\right\Vert_2 \leq 1\}$
    - prove:
      - by the definition
      - by: $P = Q^TQ$ (Affine functions preserve convexity)
  - norm ball: $\{x| \left\Vert x - x_c\right\Vert \leq r\}$
  - norm cone: $\{(x, t)| \left\Vert x \right\Vert \leq t\}$
  - polyhedron: $\{x| a^T_j x\leq b_j, c_j^T x= d_j\}$
  - simplex: a convex hull of $k + 1$ *affinely independent* points
    - we can describe the simplex as a polyhedron
  - positive semidefinite cone
    - prove:
      - by the definition
      - by: $\mathbb{S}^n_{+} = \bigcap\limits_{a\neq 0}\{x|a^Txa\succeq 0\}$

## Operations that preserve convexity

- Intersection
- Affine function: $f(x) = Ax + b, f: \mathbb{R}^n\rightarrow \mathbb{R}^m$
- Perspective function: $P(z, t) = z/t, P: \mathbb{R}^n \times \mathbb{R}_{++} \rightarrow \mathbb{R}^n$
- Linear-fractional function: $f(x) = (Ax + b)/(c^T x + d), f: \mathbb{R}^n\rightarrow \mathbb{R}^m$

We can think of affine($c=0, d>0$) and perspective functions($c = e_n, d = 0$) as special cases of linear-fractional functions.

The **image** and **inverse image** of a convex set under these three functions is convex. That is $f(S) = \{f(x)|x \in S\}$ and $f^{-1}(S) = \{x|f(x) \in S\}$ is convex if $S$ is convex.

- Prove:
  - For perspective function:
    $$
    f(\theta (x_1, t_1) + (1 - \theta)(x_2, t_2)) = \frac{\theta x_1 + (1 - \theta) x_2}{\theta t_1 + (1 - \theta) t_2} \\
    = \mu f(x_1, t_1) + (1 - \mu)f(x_2, t_2)
    $$
  - For linear-fractional function:
    Suppose $P(x)$ is the perspective function, then
    $$
    f(x) = P(Q(x, 1)^T), Q = \begin{pmatrix}
      A & b \\ c^T & d
    \end{pmatrix}
    $$

- Examples

  - Scaling: $\{\alpha x | x \in S\}$
  - Translation: $\{x + \alpha | x \in S\}$
  - Projection: $\{x_1 \in \mathbb{R}^m| (x_1, x_2)\in S, x_2 \in \mathbb{R}^n\}$
  - Cartesian product: $S_1\times S_2 = \{(x_1, x_2)|x_1\in S_1, x_2 \in S_2\}$
    - because sum of two sets(i.e. $S_1 + S_2 = \{x + y|x\in S_1, y \in S_2\}$) is convex
  - Solution set of linear matrix inequality: $\{x|A(x) = x_1A_1 + \cdots + x_nA_n \succeq B\}$
    - the inverse image of the positive semidefinite cone under the affine function: $f(x) = A(x) - B$
  - Hyperbolic cone: $\{x| x^TPx \leq (c^Tx)^2, c^Tx > 0\}, P\in \mathbb{S}^n_{+}$
    - the inverse image of the second-order cone under the affine function: $f(x) = (P^\frac{1}{2}x, c^Tx)$

## Separating and supporting hyperplanes

**Separating hyperplane theorem**: Suppose $C$ and $D$ are convex and $C \cap D = \emptyset$, then there exist $a \neq 0$ and $b$ such that $\forall x \in C, a^Tx \leq b$ and $\forall x \in D, a^Tx \geq b$. The hyperplane $\{x|a^Tx = b\}$ is called *separating hyperplane* and it separates the two sets.

- **Example**: If $C$ is convex and $D$ is affine(i.e. $D = \{Fu + g| u \in \mathbb{R}^m\}$) and $C \cap D = \emptyset$, then there exists a separating hyperplane $\{x|a^Tx = b\}$. So $a^T(Fu + g) \leq b, \forall u \in \mathbb{R}^m$ or $a^T(Fu + g) \geq b, \forall u \in \mathbb{R}^m$. Since a linear function is unbounded, $a^TF = 0$.
- Stronger condition: Suppose $C$ and $D$ are convex and $C \cap D = \emptyset$, if there exist $a \neq 0$ and $b$ such that $\forall x \in C, a^Tx < b$ and $\forall x \in D, a^Tx > b$, we call the hyperplane strictly separates the two sets.
- The **converse** of this theorem is not true. But we can add some conditions(e.g. one of the set is **open**).

**Supporting hyperplanes**: Suppose $x_0$ is a point in the boundary of a set $C$, if $\forall x\in C, a^Tx \leq a^Tx_0, a \neq 0$, then the hyperplane $\{x|a^Tx = a^Tx_0\}$ is called a *supporting hyperplanes* to $C$ at $x_0$.

**Supporting hyperplane theorem**: For any convex set $C\neq \emptyset$ and any $x_0\in \mathbf{bd}C$, there exists a supporting hyperplane to $C$ at $x_0$.
- The **converse** of this theorem is also not true. It is true when the set is closed and has nonempty interior.