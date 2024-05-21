# Graph

Notation: `G(V, E)`, `G(N, E)`, `G(N, A)`

Classification:

- Undirected and Directed Graph

- Sparse and Dense Graph

Representation:

- Adjacency Matrix (for Dense Graph)

- Adjacency Linked List (for Sparse Graph)

## Breadth First Search

`G(V, E)`, sourse `s`

BFS explores edges of `G` to discover every vertex that is reachable from `s` and provide the **shortest path from `s`** to every reachable vertex.

Notation:
- O: not discovered
- ◻︎: all neighbours having been explored
- △: in between

For each iteration, we update the node's previous node and distance from the source. (`u.d`: distance from s to u; `u.π`: privious node)

If `(u, v) in E` and `u` is a ◻︎, `v` is either △ or ◻︎; if `(u, v) in E` and `u` is a △, `v` can be △ or ◻︎ or O;

```
BFS(G, s)
	for each vertex in G.V - {s}
		u.shape = O
		u.d = inf
		u.π = NIL
	s.shape = △
	s.d = 0
	s.π = NIL
	Q = {}
	ENQUEUE(Q, s)
	while IS-EMPTY(Q) != 0
		u = DEQUEUE(Q)
		for each v in u.Adj
			if v.shape == O
				v.shape = △
        v.d = u.d + 1
        v.π = u
        ENQUEUE(Q, v)
    u.shape = ◻︎
```

Running time complexity:
- Line 2-5: `Theta(|V|)`
- Line 6-10: `Theta(1)`
- Line 11-19: `O(2|E|)`
So, the total time is `O(|V| + |E|)`.

```
PRINT-PATH(G, s, v)
	if v == s
		print s
	else if v.π == NIL
		print "No Path!"
	else
		PRINT-PATH(G, s, v.π)
    print v
```

Define the predecessor subgraph `G_π = (V_π, E_π)` where `V_π = {v in V: v.π ≠ NIL or v = s` and `E_π = {(v.π,v): v in V_π - {s}}`.

## Depth First Search

DFS searches deeper whenever possible.

1. Search unexplored edges of most recently discovered vertex.
2. When all `v`'s edges have been explored, the search "back track" to explore edges leaving vertex from which `v` was discovered.

### Two time stamp:

`v.d` records when `v` is first discovered; `v.f` records when `v` is finished.

```
DFS(G)
  for each vertex u in G.V
    u.shape = O
    u.π = NIL
  time = 0
  for each vertex u in G.V
    if u.shape == O
      DFS-VISIT(u)
		
DFS-VISIT(u)
	u.shape = △
	time ++
	u.d = time
	for each v in u.Adj
		if v.shape == O
			v.π = u
			DFS-VISIT(v)
	u.shape = ◻︎
  time = time + 1
	u.f = time
  
Topological-SORT(G)
	call DFS(G) to compute finish times u.f for each vertex u
  As each vertex is finished, insert it into the front of linked list
  return the linked list
```

Back edge: connecting `u` to an ancestor `v` in a DFS tree. Acyclic graph is a graph without back edges.

Strongly connected component of a directed graph `G` is a maximal set of verticles `C` (belongs to `V`) such that for each pair of vertices `u` and `v` in `C`, we have both `u->v` and `v->u`.

Transpose of a graph G, denoted by `G^T = (V, E^T)`, where `E^T = {(u, v): (v, u) in G}`.

```
SCC(G)
	call DFS(G) to compute u.f for each u in U
	compute G^T
	call DFS(G^T) but in the main loop of DFS consider the vertices in order of decreasing u.f
	output the vertices in the DFS forest formed in Line4
```

## Minimum Spanning Tree

- Shortest path problem –– polynomial algorithm
- Broadcast problem –– polynomial algorithm
- Multi-cast problem –– NP hard

### Greedy Algorithm

Given an undirected graph `G(V, E)`, for each edge there is a cost `w(u, v)`. Find a tree that connecting all the vertices with minimum cost.

```
GENERIC-SPANNING-TREE(G)
	A = ø
	while A does not form a spanning tree
		find an edge (u, v) that is safe for A
		A = A U {(u, v)}
	return A
```

```
KRUSTKAL'S ALGORITHM(G, W)
	A = ø
	for each vertex v in V
		makeset(v) // create a set with only one element v
	sort edges into non-decreasing order by W
	for each edge (u, v) in E, taken in non-decreasing order
		if FIND-SET(u) ≠ FIND-SET(v)
			A = A U {(u, v)}
			union(u, v)
	return A
```

```
PRIM'S ALGORITHM(G, W, r)
	for each u in V
		u.key = inf // v.key is the min weight of any edge connecting v to any vertex in a tree
		u.π = NIL
	r.key = 0
	Q = V
	while Q ≠ ø
		u = EXTRACT-MIN(Q)
		for each v in u.Adj
			if v in Q and w(u, v) < v.key
				v.π = u
				v.key = w(u, v)
```