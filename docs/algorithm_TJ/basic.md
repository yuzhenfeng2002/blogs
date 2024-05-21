# Basic Data Structures

## Stack

Last in, first out. (LIFO)

`TOP[S]` is the **index** of the most recently inserted element. When `TOP[S] = 0`, the stack is empty.

```
IS-EMPTY(S)
    if TOP[S] == 0
        return true
    else
        return false

PUSH(S, x)
    TOP[S] = TOP[S] + 1
    S[TOP[S]] = x

POP(S)
    if IS-EMPTY(S)
        error "Stack Underflow"
    else
        TOP[S] = TOP[S] - 1
        return S[TOP[S] + 1]
```

## Queue

First in, first out. (FIFO)

A queue has a `HEAD` and a `TAIL`. Insertion is called enqueue and deletion is called dequeue.

```
ENQUEUE(Q, x)
    Q[TAIL(Q)] = x
    if TAIL(Q) = Q.length
        TAIL(Q) = 1
    else
        TAIL(Q) = TAIL(Q) + 1

DEQUEUE(Q)
   x = Q[HEAD(Q)]
   if HEAD(Q) == Q.length
      HEAD(Q) = 1
   else
      HEAD(Q)  = HEAD(Q) + 1
   return x
```

# Hash Table

A hash table is like a dictionary. Under certain assumption, the searching operation is `O(1)`.

## Direct-Address Table

An array of objects with the same type is denoted by `T[0,\cdots,m-1]`. If we want to store a set of objects with a unique `k`, we will put element of key `k` into slot `k` of the array.

If the set contains no element with key `k`, the `T[k] = NIL`.

### Drawbacks

Waste memory when `\max{K} >> |K|`. When `\max{K}` is large, allocating memory is expensive.

## Hash Table

In Hash Table, each element is stored in slot `h(k)`. `h(k)` is called a hash function, computed from k. h function maps a universe of K into slot of hash table `T[0,\cdots,m-1]`.

When `h(k_1)=h(k_2)`, two keys have the same slot which is called a collision, which cannot be avoided.

### Solution

By Chaining:

```
CHAINED-HASH-INSERT(T, x)
    insert x at the head of list T[h(x.key)]
 
CHAINED-HASH-SEARCH(T, k)
    search for an element with k in list T[x.key]

CHAINED-HASH-DELETE(T, x)
	delete from the list T[h(x.key)]
```

The expected length of each list: `E[n_j] = \frac{n}{m} = \alpha`. So, the time complexity of the three function are `O(1), O(1 + \alpha), O(1 + \alpha)`.

## Hash Function

### Hash by Division

```latex
h(k)={k} \mod {m}
```

It is better to select a primer as `m`.

### Hash by Multiplication

```latex
h(k)=\lfloor m(kA \mod 1)\rfloor
```

A suitable `A` can be `\frac{\sqrt{5} - 1}{2}`.

An appliction of Hash Table is **Sim Hash**. It turns **key words** on a web page to a specific vector of 0-1.

## Binary Search Tree

### Binary search tree property

Let x be a node in a binary search tree. If y is a node in the left subtree of x, then y.key <= x.key. If y is a node in the right subtree of x, then y.key >= x.key.

```
TREE-WALK(T, x)
    IF x != NIL
        TREE-WALK(T, x.left)
        print x.key
        TREE-WALK(T, x.right)
```

It takes `Theta(n)` time.

### Dynamic set operations

```
TREE-MINIMUM(x)
    while x.left != NIL
        x = x.left
    return x

ITERATIVE-TREE-SEARCH(x, k)
    while x != NIL and k != x.key
        if k < x.key
            x = x.left
        else
            x = x.right
    return x

TREE-SUCCESSOR(x)
    if x.right == NIL
        return TREE-MINIMUM(x.right)
    y = x.p
    while y != NIL and x == y.right
        x = y
        y = y.p
    return y
```

Each one runs in O(h) time on a binary search tree of height h.

### Insertion and Deletion

```
TREE-INSERT(T, z)
    y = NIL
    x = T.root
    while x != NIL
        y = x
        if  z.key < x.key
            x = x.left
        else
            x = x.right
    z.p = y
    if y == NIL
        T.root =  z
    elseif z.key < y.key
        y.left = z
    else
        y.right = z
```

There are 3 cases for tree deletion:

1. If `z` has no child, just remove it.
2. If `z` has one child, we splice out `z`.
3. If `z` has two child, we splice out its successor `y`, which has at most one child, then replace `z.key` with `y.key`.

```
TREE-DELETE(T, z)
	if z.left == NIL or z.right == NIL
		y = z
	else
		y = TREE-SUCCESSOR(z)
    // y: node to replace z

	if y.left != NULL
		x = y.left
	else
		x = y.right
	if x != NULL
		x.p = y.p
	if y.p == NULL
		T.root = x
	else if y == y.p.left
		y.p.left = x
	else
		y.p.right = x
    // delete node y

	if y != z
		z.key = y.key
```

Each one runs in O(h) time on a binary search tree of height h too.