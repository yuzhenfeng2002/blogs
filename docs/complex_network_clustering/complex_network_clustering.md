# Complex Network Clustering

## Spectral clustering

### Preliminaries
The *input* of spectral clustering is a **Similarity Matrix** $S\in \mathbb{S}^n$ and number of clusters $k$. Every element of the similarity matrix $s_{ij}$ measures the pairwise similarity between node $i$ and $j$ on a graph $G(V, E)$. Suppose the weight of the edge between $i$ and $j$ is $w_{ij}$, then $s_{ij} = S(w_{ij})$.

Intuitively, we can minmize the cut between nodes inside and outside a cluster. That is:
$$
\min \frac{1}{2}\sum\limits_{i=1}^k W(A_i, \overline{A_i})
$$
where $A_i$ denotes a set of nodes in the $i$th cluster, $\overline{A_i}$ denotes those are not in the cluster and:
$$
W(A, B) = \sum\limits_{i\in A, j\in B}w_{ij}
$$
When $k = 2$, it is a mincut problem. In practice the solution sometimes simply separates one individual node from the others. So we adjust our *objective* as minimizing **RatioCut** or **Ncut**:
$$
\mathrm{RatioCut}(A_1, \cdots, A_k) = \frac{1}{2}\sum\limits_{i=1}^k \frac{W(A_i, \overline{A_i})}{|A_i|}\\
\mathrm{Ncut}(A_1, \cdots, A_k) = \frac{1}{2}\sum\limits_{i=1}^k \frac{W(A_i, \overline{A_i})}{\mathrm{vol}(A_i)}
$$
where $|A_i|$ means number of nodes in the $i$th cluser and $\mathrm{vol}(A_i) = \sum\limits_{k \in A_i}\sum\limits_{j = 1}^n w_{kj}$.

### Minimize RatioCut approximately
Define $H\in \mathbb{R}^{n\times k}$ as an indicator matrix which indicates whether node $i$ belongs to cluster $j$. That is:
$$
h_{ij} = \begin{cases}
    \frac{1}{\sqrt{|A_j|}}, v_i \in A_j\\
    0, \mathrm{otherwise}
\end{cases}
$$
Then we can get that:
$$
\mathrm{RatioCut}(A_1, \cdots, A_k) = \mathbf{Tr}(H^TLH)
$$
where $L$ is the **Laplacian Matrix** of the graph:
$$
L = D - W = \mathbf{diag}(\sum\limits_{j=1}^n w_{1j}, \cdots, \sum\limits_{j=1}^n w_{nj}) - (w_{ij})_{n\times n}
$$
So the optimization problem is:
$$
\min \mathbf{Tr}(H^TLH)\\
\mathrm{s.t.\ } H = (h_{ij})_{n\times k}, h_{ij} = \begin{cases}
    \frac{1}{\sqrt{|A_j|}}, v_i \in A_j\\
    0, \mathrm{otherwise}
\end{cases}
$$
This integer programming problem is a NP hard problem. But we can relax it by allowing $h_{ij}$ to take any real values. Then the relaxed problem becomes:
$$
\min \mathbf{Tr}(H^TLH)\\
\mathrm{s.t.\ } H^TH = I
$$
where $H^TH = I$ can be checked under our definition of $H$.

It is a standard **Trace Minimization Problem** and according to **Rayleigh-Ritz Theorem** the solution is given by choosing the first $k$ eigenvectors of $L$ to form $H$.

### Minimize Ncut approximately

Similar to the techniques used for minimize RatioCut, define:
$$
h_{ij} = \begin{cases}
    \frac{1}{\mathrm{vol}(A_j)}, v_i \in A_j\\
    0, \mathrm{otherwise}
\end{cases}, H = (h_{ij})_{n\times k}
$$
Then $\mathrm{Ncut}(A_1, \cdots, A_k) = \mathbf{Tr}(H^TLH)$. So the optimization problem is:
$$
\min \mathbf{Tr}(H^TLH)\\
\mathrm{s.t.\ } H = (h_{ij})_{n\times k}, h_{ij} = \begin{cases}
    \frac{1}{\mathrm{vol}(A_j)}, v_i \in A_j\\
    0, \mathrm{otherwise}
\end{cases}
$$
Relax the problem by allowing $h_{ij}$ to take any real values. Then the relaxed problem becomes:
$$
\min \mathbf{Tr}(H^TLH)\\
\mathrm{s.t.\ } H^TDH = I
$$
Let $T = D^{1/2}H$, then the problem becomes:
$$
\min \mathbf{Tr}(T^TD^{-1/2}LD^{-1/2}T)\\
\mathrm{s.t.\ } T^TT = I
$$

It is also a standard Trace Minimization Problem and the solution is given by choosing the first $k$ eigenvectors of $L_{\mathrm{sym}} = D^{-1/2}LD^{-1/2}$ as cloumns of $T$. Then we get $H = D^{-1/2}T$.

Finally, we should re-convert the real values of $H$ into discrete partition. Here we can use k-means algorithms on the rows of $H$.

### Algorithms

Then we get our algorithms for complex network clustering:

```
Unnormalized-Spectral-Clustering(S, k)
    Compute the Laplacian Matrix L
    Let H contain the first k eigenvectors of L as columns
    labels = K-Means(rows of H)
    return labels
```

```
Normalized-Spectral-Clustering(S, k)
    Compute the Normalized Laplacian Matrix Lsym
    Let T contain the first k eigenvectors of Lsym as columns
    H = D^(-1/2)T
    labels = K-Means(rows of H)
    return labels
```

## References

- von Luxburg, U. A tutorial on spectral clustering. Stat Comput 17, 395–416 (2007).
  - https://doi.org/10.1007/s11222-007-9033-z