# Running Time Analysis

## Insertion sort

Pesudocode:
```
INSERTION-SORT(A)
for j from 2 to A.length
    key <- A[j]
    i <- j - 1
    while i > 0 and A[i] > key
        A[i + 1] <- A[i]
        i <- i - 1
    A[i + 1] <- key
```

## Analysis of algorithms

Analysis of algorithms is to find the relation between input and its running time.

running time: the number of primitive operations/steps

> for `INSERTION-SORT`:
>
> Worst-case analysis: `1+2+\cdots+(n-1) = \frac{n(n-1)}{2} = \frac{n^2}{2} - \frac{n}{2}`;
>
> Order of growth: `\Theta(n^2)`
>
> Average case analysis: `\frac{n(n-1)}{4} = \frac{n^2}{4} - \frac{n}{4}`
>
> Order of growth: `\Theta(n^2)`

Formally,

`\Theta(g(n)) = f(n)`: there exists positive constant `c_1`, `c_2` and `n_0` such that `0 \leq c_1g(n) \leq f(n) \leq c_2g(n)` for all `n \geq n_0`;

`O(g(n)) = f(n)`: there exists positive constant `c` and `n_0` such that `0 \leq f(n) \leq cg(n)` for all `n \geq n_0`.

If we want to prove `\frac{n^2}{2}-3n=\Theta(n^2)`, we can assume `c_1n^2 \leq \frac{n^2}{2}-3n \leq c_2n^2`.

```latex
c_1n^2 \leq \frac{n^2}{2}-3n \leq c_2n^2 \\
\Rightarrow c_1 \leq \frac{1}{2}-\frac{3}{n} \leq c_2
```

Let `c_1 = \frac{1}{4}, c_2 = \frac{1}{2}`.

# Recursive and Recurrences

## Divide and conquer

Many Algorithms are recursive in structure. These algorithm call itself one or more times to deal with the closely related **subproblem** with similar structure.

Steps:
- Divide: Break the problem into several subproblems that are similar to the original problem, but smaller in size.
- Conquer: Solve the subproblem recursively. when subproblem is small enough, solve the subproblem straightforward.
- Combine: Combine the solutions to the subproblems to the original problem.

## Merge sort

- Divide: Divide the n-element sequence into 2 sub sequences of n/2 elements each.
- Conquer: Sort 2 sub sequences using merge sort.
- Combine: Merge the **sorted** sub sequences to provide the answer.

Pesudocode:
```
MERGE(A, p, q, r)
    n1 <- q - p + 1
    n2 <- r - q
    for i <- 1 to n1
        L[i] <- A[p + i - 1]
    for j <- 1 to n2
        R[i] <- A[q + j]
    L[n1 + 1] <- infinity
    R[n2 + 1] <- infinity
    i <- 1, j <- 1
    for k from p to r
        if L[i] <= R[j]
            A[k] <- L[i]
            i <- i + 1
        else 
            A[k] <- R[j]
            j <- j + 1

MERGESORT(A, p, r)
    if p < r
        q <- FLOOR(((p + r) / 2)
        MERGESORT(A, p, q)
        MERGESORT(A, q + 1, r)
        MERGE(A, p, q, r)
```

## Runtime of Merge Sort

```latex
T(n) = \begin{cases}
    1,n=1\\
    2T(\frac{n}{2}) + \Theta(n), n\geq 2
\end{cases}
```

### Induction/ Substitution Method
- Step 1: 
  When `n_0 = 2`, there exists a constant `c` such that `T(n_0) \leq c(2\lg2)`.
- Step 2:
  ```latex
  T(n) = 2T(\frac{n}{2}) + \Theta(n)\\
  \leq 2(c\frac{n}{2}\lg \frac{n}{2}) + cn\\
  = cn\lg n - cn\lg2 + cn\\
  = cn\lg n.
  ```

## Master Theorem

Let `a \geq 1`, `b > 1` be constants, `f(n)` be a function and `T(n)` be defined by recurrence:
```T(n) = aT(\frac{n}{b}) + f(n)```

1. If `f(n) = O(n^{\log _b a - \epsilon})` for some constant `\epsilon > 0`, then `T(n) = \Theta(n^{\log _b a})`;
2. If `f(n) = \Theta (n^{\log _b a})`, then `T(n) = \Theta(n^{\log _b a}\lg n)`;
3. If `f(n) = \Omega (n^{\log _b a + \epsilon})` for some constant `\epsilon > 0`, and if `af(\frac{n}{b})\leq cf(n)` for some constant `c < 1` and sufficiently large `n`, then `T(n) = \Theta(f(n))`.

# Heap and Heapsort

## Heap

A heap is an array with special structure.

The root of a heap is A[1]. For any element i in the array, it **may** have a parent, a left child and a right child.

```
PARENT(i)
    return floor(i/2)
LEFT(i)
    return 2i
RIGHT(i)
    return 2i+1
```

Heap property:

- Max Heap: `A[PARENT(i)]` >= `A[i]`
- Min Heap: `A[PARENT(i)]` <= `A[i]`

## Build a heap

`MAX-HEAPIFY(A, i)` is based on the assumption that 2 trees rooted at `LEFT(i)` and `RIGHT(i)` are heaps.

```
MAX-HEAPIFY(A, i)
    l = LEFT(i); r = RIGHT(i)
    IF l <= A.size() and A[l] > A[i]
        largest = l
    ELSE largest = i
    IF r <= A.size() and A[r] > A[i]
        largest = r
    IF largest != i
        exchange A[largest] with A[i]
        MAX-HEAPIFY(A, largest)
```

Runtime: `O(\lg n)`

```
BUILD-HEAP(A)
    FOR i from floor(A.leangth / 2) down to 1
        MAX-HEAPIFY(A, i)
```

Runtime: `O(n) = O(n\lg n)`

Hint (refer to *CLRS* for details):

```latex
T(n) = \sum \limits _{h=0}^{\lfloor\lg n\rfloor} \lfloor\frac{n}{2^{h+1}}\rfloor O(h)
```

## Heap sort

```
HEAP-SORT(A)
    BUILD-HEAP(A)
    FOR i from A.length down to 2
        EXCHANGE A[1] with A[i]
        A.size = A.size - 1
        MAX-HEAPIFY(A, 1)
```

Runtime: `O(n + n\lg n) = O(n\lg n)`

# Quick Sort

1. Divide: Divide A[p, ..., r] into 2 subs: A[p, ..., q-1], A[q+1, ..., r] s.t. elements in the formal smaller than A[q] and  the latter larger than A[q].
2. Conquer: Sort 2 subs recursively by calling `Quicksort`.
3. Combine: Do nothing.

```
QUICKSORT(A, p, r)
	if p < r
		q = PARTITION(A, p, r)
		QUICKSORT(A, p, q - 1)
		QUICKSORT(A, q + 1, r)

PARTITION(A, p, r)
	x = A[r]
	i = p - 1
	for j = p to r - 1
		if A[j] <= x
			i = i + 1
			exchange A[i] with A[j]
	exchange A[i + 1] with A[r]
	return i + 1
```

## Time Complexity

```latex
T(n) = T(q) + T(n-q-1) + \Theta(n)
```

### Worst case

```latex
T(n) = T(n - 1)  + \Theta(n)\\
T(n) = \Theta(n^2)
```

### Best case

```latex
T(n) = 2T(\frac{n}{2}) + \Theta(n)\\
T(n) = \Theta(n\log n)
```
