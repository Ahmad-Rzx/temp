"""
graphs.py  --  BFS/DFS, shortest paths, topo sort, SCC, 2-SAT, bridges, MST, trees/LCA, max flow
===============================================================================
CONVENTION: every algorithm uses 0-indexed nodes. Read 1-indexed input with build_graph(...,
one_indexed=True) (it subtracts 1) and add 1 again when printing.
adj[u] = [v, ...]   (unweighted)    or    adj[u] = [(v, w), ...]   (weighted)
All DFS-like routines are ITERATIVE (no recursion-limit problems).
"""
from collections import deque
from heapq import heappush, heappop

INF = float("inf")


def build_graph(n, edges, directed=False, weighted=False, one_indexed=True):
    """edges: [(u, v)] or [(u, v, w)] -> adjacency list of size n."""
    off = 1 if one_indexed else 0
    adj = [[] for _ in range(n)]
    for e in edges:
        u, v = e[0] - off, e[1] - off
        if weighted:
            w = e[2]
            adj[u].append((v, w))
            if not directed:
                adj[v].append((u, w))
        else:
            adj[u].append(v)
            if not directed:
                adj[v].append(u)
    return adj


# ═══════════════════════════════ BFS FAMILY ═══════════════════════════════════
# WHEN:  shortest path in an UNWEIGHTED graph / grid, "minimum number of moves/steps",
#        level order, connected components, multi-source spreading (fire, rot, flood).
# MODIFY: pass a list of sources for multi-source BFS. For state-space BFS (e.g. position + keys
#         held, or (x, y, direction)) use tuples/encoded ints as nodes and a dict/array for dist.
def bfs(adj, src):
    """Returns (dist, parent). src may be an int or a list of sources. dist = -1 if unreachable."""
    n = len(adj)
    dist, parent = [-1] * n, [-1] * n
    srcs = [src] if isinstance(src, int) else list(src)
    q = deque(srcs)
    for s in srcs:
        dist[s] = 0
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                parent[v] = u
                q.append(v)
    return dist, parent


def restore_path(parent, t):
    """Path from the BFS/Dijkstra source to t (list of nodes), using the parent array."""
    path = []
    while t != -1:
        path.append(t)
        t = parent[t]
    return path[::-1]


def components(adj):
    """Returns (comp_id list, number_of_components) for an undirected graph."""
    n = len(adj)
    comp = [-1] * n
    c = 0
    for s in range(n):
        if comp[s] >= 0:
            continue
        comp[s] = c
        st = [s]
        while st:
            u = st.pop()
            for v in adj[u]:
                if comp[v] < 0:
                    comp[v] = c
                    st.append(v)
        c += 1
    return comp, c


def grid_bfs(grid, sr, sc, wall="#"):
    """Shortest steps from (sr, sc) in a grid of strings/lists, 4-directional. -1 = unreachable.
    MODIFY: add 4 diagonal moves for 8-directional; change the `wall` test for other blockers."""
    R, C = len(grid), len(grid[0])
    dist = [[-1] * C for _ in range(R)]
    dist[sr][sc] = 0
    q = deque([(sr, sc)])
    while q:
        r, c = q.popleft()
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            nr, nc = r + dr, c + dc
            if 0 <= nr < R and 0 <= nc < C and grid[nr][nc] != wall and dist[nr][nc] < 0:
                dist[nr][nc] = dist[r][c] + 1
                q.append((nr, nc))
    return dist


def bfs01(adj, src):
    """0-1 BFS: weights are only 0 or 1. adj[u] = [(v, w)]. O(V+E). WHEN: 'cost 1 to flip/turn, 0 to continue'."""
    dist = [INF] * len(adj)
    dist[src] = 0
    dq = deque([src])
    while dq:
        u = dq.popleft()
        for v, w in adj[u]:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                if w == 0:
                    dq.appendleft(v)
                else:
                    dq.append(v)
    return dist


def is_bipartite(adj):
    """Returns (True, colors) or (False, None). WHEN: 'split into two groups with no edge inside', odd cycle."""
    n = len(adj)
    col = [-1] * n
    for s in range(n):
        if col[s] >= 0:
            continue
        col[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if col[v] < 0:
                    col[v] = col[u] ^ 1
                    q.append(v)
                elif col[v] == col[u]:
                    return False, None
    return True, col


# ═════════════════════════════ SHORTEST PATHS ═════════════════════════════════
# DIJKSTRA  WHEN: non-negative weights, single source.  TIME O((V+E) log V).
# MODIFY: state-space Dijkstra = same code with tuple nodes; to count shortest paths also keep
#         ways[v] (add ways[u] when dist ties, reset when strictly better).
def dijkstra(adj, src):
    """adj[u] = [(v, w)]. Returns (dist, parent)."""
    n = len(adj)
    dist = [INF] * n
    parent = [-1] * n
    dist[src] = 0
    pq = [(0, src)]
    while pq:
        d, u = heappop(pq)
        if d > dist[u]:
            continue                             # stale entry
        for v, w in adj[u]:
            nd = d + w
            if nd < dist[v]:
                dist[v] = nd
                parent[v] = u
                heappush(pq, (nd, v))
    return dist, parent


def bellman_ford(n, edges, src):
    """edges: [(u, v, w)] directed. Handles NEGATIVE weights. Returns (dist, has_negative_cycle).
    WHEN: negative edges, detect negative cycle (arbitrage). O(VE)."""
    dist = [INF] * n
    dist[src] = 0
    for i in range(n):
        changed = False
        for u, v, w in edges:
            if dist[u] < INF and dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                changed = True
        if not changed:
            return dist, False
    return dist, True


def floyd_warshall(n, edges, directed=True):
    """All-pairs shortest paths. WHEN: n <= ~400 (O(n^3)). Returns matrix d (INF if unreachable).
    Negative cycle exists iff some d[i][i] < 0."""
    d = [[INF] * n for _ in range(n)]
    for i in range(n):
        d[i][i] = 0
    for u, v, w in edges:
        if w < d[u][v]:
            d[u][v] = w
        if not directed and w < d[v][u]:
            d[v][u] = w
    for k in range(n):
        dk = d[k]
        for i in range(n):
            dik = d[i][k]
            if dik == INF:
                continue
            d[i] = [x if x <= dik + y else dik + y for x, y in zip(d[i], dk)]
    return d


# ═════════════════════════ TOPOLOGICAL ORDER / SCC / 2-SAT ════════════════════
# WHEN:  dependencies ("A must come before B"), detect cycle in directed graph, DP on a DAG.
# MODIFY: use a heap instead of a deque to get the LEXICOGRAPHICALLY SMALLEST order.
def topo_sort(adj):
    """Kahn's algorithm. Returns an order, or None if the graph has a cycle."""
    n = len(adj)
    indeg = [0] * n
    for u in range(n):
        for v in adj[u]:
            indeg[v] += 1
    q = deque(i for i in range(n) if indeg[i] == 0)
    order = []
    while q:
        u = q.popleft()
        order.append(u)
        for v in adj[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                q.append(v)
    return order if len(order) == n else None


def topo_sort_lex_smallest(adj):
    n = len(adj)
    indeg = [0] * n
    for u in range(n):
        for v in adj[u]:
            indeg[v] += 1
    h = [i for i in range(n) if indeg[i] == 0]
    order = []
    while h:
        u = heappop(h)
        order.append(u)
        for v in adj[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                heappush(h, v)
    return order if len(order) == n else None


def scc(adj):
    """Strongly connected components (Kosaraju, iterative). Returns (comp, count).
    comp ids are in TOPOLOGICAL order of the condensation DAG (comp 0 has no incoming edges
    from other comps). WHEN: 'mutually reachable groups', make graph strongly connected, DP over condensation."""
    n = len(adj)
    radj = [[] for _ in range(n)]
    for u in range(n):
        for v in adj[u]:
            radj[v].append(u)
    seen = [False] * n
    order = []
    for s in range(n):
        if seen[s]:
            continue
        seen[s] = True
        stack = [(s, 0)]
        while stack:
            u, i = stack.pop()
            if i < len(adj[u]):
                stack.append((u, i + 1))
                v = adj[u][i]
                if not seen[v]:
                    seen[v] = True
                    stack.append((v, 0))
            else:
                order.append(u)
    comp = [-1] * n
    c = 0
    for s in reversed(order):
        if comp[s] >= 0:
            continue
        comp[s] = c
        st = [s]
        while st:
            u = st.pop()
            for v in radj[u]:
                if comp[v] < 0:
                    comp[v] = c
                    st.append(v)
        c += 1
    return comp, c


class TwoSat:
    """2-SAT. Variables 0..n-1. add_clause(i, a, j, b) means (x_i == a) OR (x_j == b), a,b booleans.
    Force x_i = a with add_clause(i, a, i, a). solve() -> list of bools, or None if unsatisfiable.
    WHEN: 'each item picks one of two options, some pairs conflict'."""
    def __init__(self, n):
        self.n = n
        self.adj = [[] for _ in range(2 * n)]

    @staticmethod
    def _id(i, val):
        return 2 * i + (0 if val else 1)

    def add_clause(self, i, a, j, b):
        A, B = self._id(i, a), self._id(j, b)
        self.adj[A ^ 1].append(B)
        self.adj[B ^ 1].append(A)

    def solve(self):
        comp, _ = scc(self.adj)
        res = []
        for i in range(self.n):
            if comp[2 * i] == comp[2 * i + 1]:
                return None
            res.append(comp[2 * i] > comp[2 * i + 1])
        return res


# ═══════════════════ BRIDGES & ARTICULATION POINTS (iterative) ════════════════
# WHEN:  "edge whose removal disconnects the graph", "critical connections", cut vertices.
# edges: [(u, v)] 0-indexed, multi-edges allowed. Returns (bridge_edge_indices, cut_vertices).
def bridges_and_cut_vertices(n, edges):
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
    tin = [-1] * n
    low = [0] * n
    pe = [-1] * n                                # parent edge id
    it = [0] * n
    timer = 0
    bridges, is_cut = [], [False] * n
    for root in range(n):
        if tin[root] != -1:
            continue
        tin[root] = low[root] = timer
        timer += 1
        stack = [root]
        root_children = 0
        while stack:
            u = stack[-1]
            if it[u] < len(adj[u]):
                v, eid = adj[u][it[u]]
                it[u] += 1
                if eid == pe[u]:
                    continue
                if tin[v] == -1:
                    pe[v] = eid
                    tin[v] = low[v] = timer
                    timer += 1
                    stack.append(v)
                else:
                    low[u] = min(low[u], tin[v])
            else:
                stack.pop()
                if stack:
                    p = stack[-1]
                    low[p] = min(low[p], low[u])
                    if low[u] > tin[p]:
                        bridges.append(pe[u])
                    if p == root:
                        root_children += 1
                    elif low[u] >= tin[p]:
                        is_cut[p] = True
        if root_children > 1:
            is_cut[root] = True
    return bridges, [i for i in range(n) if is_cut[i]]


# ═══════════════════════ MINIMUM SPANNING TREE ════════════════════════════════
# WHEN:  connect all nodes with min total weight; "minimize the maximum edge" (MST minimizes it too).
def kruskal(n, edges):
    """edges: [(w, u, v)] 0-indexed. Returns (total_weight, used_edges). Graph disconnected -> forest."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    total, used = 0, []
    for w, u, v in sorted(edges):
        a, b = find(u), find(v)
        if a != b:
            parent[a] = b
            total += w
            used.append((w, u, v))
    return total, used


def prim(adj):
    """adj[u] = [(v, w)], connected graph. Returns total MST weight. Good for dense graphs."""
    n = len(adj)
    seen = [False] * n
    pq = [(0, 0)]
    total = 0
    while pq:
        w, u = heappop(pq)
        if seen[u]:
            continue
        seen[u] = True
        total += w
        for v, ww in adj[u]:
            if not seen[v]:
                heappush(pq, (ww, v))
    return total


# ═════════════════════════════════ TREES ══════════════════════════════════════
# WHEN:  anything on trees: subtree sums/sizes, depth, diameter, distances, ancestors, LCA.
# PATTERN: get BFS order + parent, then do bottom-up DP with `for u in reversed(order)`:
#     order, parent, depth = tree_order(adj)
#     size = [1]*n
#     for u in reversed(order):
#         if parent[u] >= 0: size[parent[u]] += size[u]
def tree_order(adj, root=0):
    """Unweighted adjacency. Returns (bfs_order, parent, depth)."""
    n = len(adj)
    parent = [-1] * n
    depth = [0] * n
    order = [root]
    for u in order:
        for v in adj[u]:
            if v != parent[u]:
                parent[v] = u
                depth[v] = depth[u] + 1
                order.append(v)
    return order, parent, depth


def subtree_sizes(adj, root=0):
    order, parent, _ = tree_order(adj, root)
    size = [1] * len(adj)
    for u in reversed(order):
        if parent[u] >= 0:
            size[parent[u]] += size[u]
    return size


def tree_diameter(adj):
    """Returns (length_in_edges, endpoint_a, endpoint_b). Two BFS passes. Weighted: use dijkstra twice."""
    d0, _ = bfs(adj, 0)
    a = max(range(len(adj)), key=lambda i: d0[i])
    d1, _ = bfs(adj, a)
    b = max(range(len(adj)), key=lambda i: d1[i])
    return d1[b], a, b


def euler_tour(adj, root=0):
    """tin/tout so that subtree(u) = order[tin[u]:tout[u]] and v in subtree(u) <=> tin[u] <= tin[v] < tout[u].
    WHEN: subtree queries become RANGE queries -> combine with BIT / segment tree."""
    n = len(adj)
    tin, tout = [0] * n, [0] * n
    parent = [-1] * n
    it = [0] * n
    order = [root]
    timer = 1
    stack = [root]
    while stack:
        u = stack[-1]
        if it[u] < len(adj[u]):
            v = adj[u][it[u]]
            it[u] += 1
            if v != parent[u]:
                parent[v] = u
                tin[v] = timer
                timer += 1
                order.append(v)
                stack.append(v)
        else:
            tout[u] = timer
            stack.pop()
    return tin, tout, order


class LCA:
    """Binary lifting. WHEN: distance between tree nodes, k-th ancestor, 'is u an ancestor of v',
    path queries (combine with prefix sums from the root). O(n log n) build, O(log n) query."""
    def __init__(self, adj, root=0):
        n = len(adj)
        order, parent, depth = tree_order(adj, root)
        self.depth = depth
        self.LOG = max(1, n.bit_length())
        up0 = [p if p >= 0 else root for p in parent]
        self.up = [up0]
        for _ in range(1, self.LOG):
            prev = self.up[-1]
            self.up.append([prev[prev[v]] for v in range(n)])

    def kth_ancestor(self, v, k):
        i = 0
        while k and v is not None:
            if k & 1:
                v = self.up[i][v]
            k >>= 1
            i += 1
        return v

    def lca(self, a, b):
        if self.depth[a] < self.depth[b]:
            a, b = b, a
        a = self.kth_ancestor(a, self.depth[a] - self.depth[b])
        if a == b:
            return a
        for k in range(self.LOG - 1, -1, -1):
            if self.up[k][a] != self.up[k][b]:
                a, b = self.up[k][a], self.up[k][b]
        return self.up[0][a]

    def dist(self, a, b):
        return self.depth[a] + self.depth[b] - 2 * self.depth[self.lca(a, b)]


# ═══════════════════ MAX FLOW (Dinic) + BIPARTITE MATCHING ════════════════════
# WHEN:  max flow / min cut, "assign tasks to workers", bipartite matching, max edge-disjoint paths,
#        "minimum number of things to remove to separate S and T".
# MODIFY: model the problem as a network: source -> left side (cap 1) -> right side (cap 1) -> sink.
#         Min cut value == max flow value. Vertex capacities: split a node into in/out with an edge of cap c.
class Dinic:
    def __init__(self, n):
        self.n = n
        self.g = [[] for _ in range(n)]          # edge = [to, cap, index_of_reverse_edge]

    def add_edge(self, u, v, cap):
        self.g[u].append([v, cap, len(self.g[v])])
        self.g[v].append([u, 0, len(self.g[u]) - 1])

    def max_flow(self, s, t):
        g, n = self.g, self.n
        flow = 0
        while True:
            level = [-1] * n
            level[s] = 0
            q = deque([s])
            while q:
                u = q.popleft()
                for v, cap, _ in g[u]:
                    if cap > 0 and level[v] < 0:
                        level[v] = level[u] + 1
                        q.append(v)
            if level[t] < 0:
                return flow
            it = [0] * n
            while True:                          # iterative blocking-flow DFS
                path = []
                u = s
                while u != t:
                    advanced = False
                    while it[u] < len(g[u]):
                        e = g[u][it[u]]
                        if e[1] > 0 and level[e[0]] == level[u] + 1:
                            path.append((u, it[u]))
                            u = e[0]
                            advanced = True
                            break
                        it[u] += 1
                    if not advanced:
                        if not path:
                            break
                        level[u] = -1            # dead end
                        u, _ = path.pop()
                        it[u] += 1
                if u != t:
                    break
                f = min(g[a][i][1] for a, i in path)
                for a, i in path:
                    e = g[a][i]
                    e[1] -= f
                    g[e[0]][e[2]][1] += f
                flow += f


def bipartite_matching(n_left, n_right, pairs):
    """pairs: [(left_i, right_j)] 0-indexed. Returns maximum matching size."""
    S, T = n_left + n_right, n_left + n_right + 1
    d = Dinic(n_left + n_right + 2)
    for i in range(n_left):
        d.add_edge(S, i, 1)
    for j in range(n_right):
        d.add_edge(n_left + j, T, 1)
    for i, j in pairs:
        d.add_edge(i, n_left + j, 1)
    return d.max_flow(S, T)
