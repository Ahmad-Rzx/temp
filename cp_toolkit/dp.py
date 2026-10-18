"""
dp.py  --  classic dynamic-programming templates
===============================================================================
How to attack ANY DP problem:  (1) define state dp[...] in words, (2) write the transition,
(3) fix the base case and iteration ORDER, (4) check answer cell, (5) shrink memory if needed.
Each template says WHEN / MODIFY / TIME.
"""
from bisect import bisect_left, bisect_right
from functools import lru_cache

INF = float("inf")


# ═══════════════════════════════ KNAPSACK FAMILY ══════════════════════════════
# WHEN:  choose items with weight/cost limit and maximize value; "can we reach exactly sum S".
# TIME:  O(n * W). If W is huge but values are small, swap roles: dp[value] = min weight.
def knapsack_01(items, W):
    """items: [(weight, value)], each usable ONCE. dp[c] = best value with capacity <= c.
    MODIFY: iterate c DOWNWARDS (else you reuse items). For 'exactly W' init dp=[-INF]*(W+1), dp[0]=0."""
    dp = [0] * (W + 1)
    for w, v in items:
        for c in range(W, w - 1, -1):
            if dp[c - w] + v > dp[c]:
                dp[c] = dp[c - w] + v
    return dp[W]


def knapsack_unbounded(items, W):
    """Unlimited copies of each item. Iterate c UPWARDS."""
    dp = [0] * (W + 1)
    for w, v in items:
        for c in range(w, W + 1):
            if dp[c - w] + v > dp[c]:
                dp[c] = dp[c - w] + v
    return dp[W]


def knapsack_bounded(items, W):
    """items: [(weight, value, count)]. Binary splitting turns it into a 0/1 knapsack."""
    flat = []
    for w, v, cnt in items:
        k = 1
        while cnt > 0:
            take = min(k, cnt)
            flat.append((w * take, v * take))
            cnt -= take
            k <<= 1
    return knapsack_01(flat, W)


def subset_sum_bitset(nums, target):
    """Can some subset sum to target? Python big-int as a bitset: ~64x faster than a bool DP.
    MODIFY: `bits` has bit s set iff sum s is reachable -> read ANY reachable sum (e.g. closest to total/2)."""
    bits = 1
    for x in nums:
        bits |= bits << x
    return (bits >> target) & 1 == 1


def coin_change_min(coins, amount):
    """Minimum coins to make amount (-1 if impossible). Unbounded."""
    dp = [0] + [INF] * amount
    for a in range(1, amount + 1):
        for c in coins:
            if c <= a and dp[a - c] + 1 < dp[a]:
                dp[a] = dp[a - c] + 1
    return dp[amount] if dp[amount] < INF else -1


def coin_change_ways(coins, amount, mod=None):
    """Number of COMBINATIONS (order doesn't matter): coins in the OUTER loop.
    For ordered sequences (permutations) put the amount loop OUTSIDE."""
    dp = [1] + [0] * amount
    for c in coins:
        for a in range(c, amount + 1):
            dp[a] += dp[a - c]
            if mod:
                dp[a] %= mod
    return dp[amount]


# ═════════════════════════════ SEQUENCE DPs ═══════════════════════════════════
# LIS  WHEN: "longest increasing subsequence", "min number of removals to make sorted", chains/nesting.
# MODIFY: strict=False for non-decreasing. For decreasing negate values. Two-key nesting problems
#         (envelopes): sort by (a asc, b DESC) then LIS on b.   TIME O(n log n).
def lis_length(a, strict=True):
    tails = []
    f = bisect_left if strict else bisect_right
    for x in a:
        i = f(tails, x)
        if i == len(tails):
            tails.append(x)
        else:
            tails[i] = x
    return len(tails)


def lis_sequence(a):
    """One strictly increasing subsequence of maximum length (the elements)."""
    tails, tails_idx, prev = [], [], [-1] * len(a)
    for i, x in enumerate(a):
        j = bisect_left(tails, x)
        if j == len(tails):
            tails.append(x)
            tails_idx.append(i)
        else:
            tails[j] = x
            tails_idx[j] = i
        prev[i] = tails_idx[j - 1] if j > 0 else -1
    res, k = [], tails_idx[-1] if tails_idx else -1
    while k != -1:
        res.append(a[k])
        k = prev[k]
    return res[::-1]


def lcs(a, b):
    """Longest common subsequence: returns (length, one LCS as a list). O(len(a)*len(b)).
    MODIFY: longest palindromic subsequence = lcs(s, s[::-1]); min insertions to make palindrome = n - that."""
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(m):
            if a[i] == b[j]:
                dp[i + 1][j + 1] = dp[i][j] + 1
            else:
                dp[i + 1][j + 1] = max(dp[i][j + 1], dp[i + 1][j])
    i, j, out = n, m, []
    while i > 0 and j > 0:
        if a[i - 1] == b[j - 1]:
            out.append(a[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1
    return dp[n][m], out[::-1]


def edit_distance(a, b):
    """Levenshtein distance. MODIFY: change the three costs (insert / delete / replace)."""
    n, m = len(a), len(b)
    prev = list(range(m + 1))
    for i in range(1, n + 1):
        cur = [i] + [0] * m
        for j in range(1, m + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return prev[m]


def max_subarray(a):
    """Kadane. Returns (best_sum, start, end_inclusive). Non-empty subarray.
    MODIFY: circular version = max(kadane, total - min_subarray); allow empty -> max(0, ...)."""
    best = cur = a[0]
    bs = cs = be = 0
    for i in range(1, len(a)):
        if cur < 0:
            cur, cs = a[i], i
        else:
            cur += a[i]
        if cur > best:
            best, bs, be = cur, cs, i
    return best, bs, be


# ═════════════════════════ GRID / COUNTING DPs ════════════════════════════════
def grid_paths(grid, mod=10**9 + 7, blocked="#"):
    """Number of right/down paths from top-left to bottom-right avoiding blocked cells.
    MODIFY: for min/max path cost replace the sum by min/max of the two predecessors."""
    R, C = len(grid), len(grid[0])
    dp = [[0] * C for _ in range(R)]
    if grid[0][0] != blocked:
        dp[0][0] = 1
    for i in range(R):
        for j in range(C):
            if grid[i][j] == blocked or (i == 0 and j == 0):
                continue
            dp[i][j] = ((dp[i - 1][j] if i else 0) + (dp[i][j - 1] if j else 0)) % mod
    return dp[R - 1][C - 1]


# ════════════════════════════ BITMASK DP ══════════════════════════════════════
# WHEN:  n <= ~20 and the state is "which subset of items is already used".
#        Assignment problems, TSP, covering, "partition into groups".
# TIME:  O(2^n * n) (TSP) -> n <= ~15 comfortably in Python, ~18 in PyPy.
def tsp(dist):
    """Min cost Hamiltonian cycle through all nodes starting/ending at 0. dist: n x n matrix.
    MODIFY: drop the final return edge for a Hamiltonian PATH; dp[mask][last] semantic stays."""
    n = len(dist)
    FULL = 1 << n
    dp = [[INF] * n for _ in range(FULL)]
    dp[1][0] = 0
    for mask in range(1, FULL):
        if not mask & 1:
            continue
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c == INF:
                continue
            for nxt in range(n):
                if mask >> nxt & 1:
                    continue
                v = c + dist[last][nxt]
                nm = mask | (1 << nxt)
                if v < dp[nm][nxt]:
                    dp[nm][nxt] = v
    return min(dp[FULL - 1][i] + dist[i][0] for i in range(n))


def sos_subset_sum(a, bits):
    """Sum-over-subsets: f[mask] = sum of a[sub] over all sub that are subsets of mask. O(bits * 2^bits).
    WHEN: 'count pairs with a & b == x', 'for every mask, aggregate over its submasks/supermasks'."""
    f = list(a)
    for b in range(bits):
        bit = 1 << b
        for mask in range(1 << bits):
            if mask & bit:
                f[mask] += f[mask ^ bit]
    return f


# ════════════════════════════ DIGIT DP ════════════════════════════════════════
# WHEN:  "count numbers in [L, R] with property of digits" (digit sum, no two adjacent equal, contains 7...)
# MODIFY: answer(L, R) = f(R) - f(L-1). Add state dims in go(...) for what the property needs
#         (e.g. previous digit, a `started` flag if leading zeros must be ignored).
def count_digit_sum(N, S):
    """Count integers x in [0, N] whose digit sum equals S."""
    if N < 0:
        return 0
    digits = list(map(int, str(N)))
    n = len(digits)

    @lru_cache(maxsize=None)
    def go(pos, rem, tight):
        if rem < 0:
            return 0
        if pos == n:
            return 1 if rem == 0 else 0
        hi = digits[pos] if tight else 9
        return sum(go(pos + 1, rem - d, tight and d == hi) for d in range(hi + 1))

    return go(0, S, True)


# ═══════════════════════════ INTERVAL DP ══════════════════════════════════════
# WHEN:  merge/split a contiguous segment: "min cost to merge adjacent piles", matrix chain,
#        burst balloons, optimal BST, palindrome partitioning.   TIME O(n^3).
def min_merge_cost(a):
    """Merge adjacent piles; cost of a merge = total size. Min total cost. (Template for interval DP.)
    MODIFY: change cost(l, r) and the combine rule."""
    n = len(a)
    if n <= 1:
        return 0
    pre = [0]
    for x in a:
        pre.append(pre[-1] + x)
    dp = [[0] * n for _ in range(n)]
    for length in range(2, n + 1):
        for l in range(n - length + 1):
            r = l + length - 1
            best = INF
            for k in range(l, r):
                v = dp[l][k] + dp[k + 1][r]
                if v < best:
                    best = v
            dp[l][r] = best + pre[r + 1] - pre[l]
    return dp[0][n - 1]


# ═════════════════════════════ TREE DP ════════════════════════════════════════
# WHEN:  choose nodes/edges in a tree under constraints. Process nodes in reversed BFS order so
#        children are done before parents (no recursion). Rerooting DP = same + a second top-down pass.
def tree_max_weight_independent_set(adj, w, root=0):
    """Max total weight of a set of nodes with no two adjacent. adj unweighted tree."""
    n = len(adj)
    parent = [-1] * n
    order = [root]
    for u in order:
        for v in adj[u]:
            if v != parent[u]:
                parent[v] = u
                order.append(v)
    take = list(w)
    skip = [0] * n
    for u in reversed(order):
        p = parent[u]
        if p >= 0:
            take[p] += skip[u]
            skip[p] += max(take[u], skip[u])
    return max(take[root], skip[root])


# ═══════════════════════ LI CHAO TREE (convex hull trick) ═════════════════════
# WHEN:  dp[i] = min over j of (m_j * x_i + b_j)  -- "lines" DP, slope-trick-like optimizations.
#        Query min of many lines at a point.  Pass ALL x values you will ever query (sorted, unique).
# MODIFY: for MAX negate slopes/intercepts and the result.
class LiChao:
    def __init__(self, xs):
        self.xs = xs
        self.n = len(xs)
        self.t = [None] * (4 * max(1, self.n))

    def add_line(self, m, b):
        node, lo, hi = 1, 0, self.n - 1
        xs, t = self.xs, self.t
        while True:
            cur = t[node]
            if cur is None:
                t[node] = (m, b)
                return
            mid = (lo + hi) // 2
            xm = xs[mid]
            if m * xm + b < cur[0] * xm + cur[1]:
                t[node], (m, b) = (m, b), cur    # keep the better-at-mid line in the node
            if lo == hi:
                return
            cur = t[node]
            if m * xs[lo] + b < cur[0] * xs[lo] + cur[1]:
                node, hi = 2 * node, mid
            elif m * xs[hi] + b < cur[0] * xs[hi] + cur[1]:
                node, lo = 2 * node + 1, mid + 1
            else:
                return

    def query(self, x):
        idx = bisect_left(self.xs, x)
        node, lo, hi = 1, 0, self.n - 1
        res = INF
        while node < len(self.t) and self.t[node] is not None:
            m, b = self.t[node]
            v = m * x + b
            if v < res:
                res = v
            if lo == hi:
                break
            mid = (lo + hi) // 2
            if idx <= mid:
                node, hi = 2 * node, mid
            else:
                node, lo = 2 * node + 1, mid + 1
        return res
