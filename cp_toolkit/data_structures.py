"""
data_structures.py  --  DSU, Fenwick, segment trees, sparse table, monotonic stack/deque,
                        sorted list, xor trie, removable heap, running median
===============================================================================
All indices are 0-based. Ranges are half-open [l, r) unless the docstring says otherwise.
Each class is self-contained: copy only what you need.
"""
from bisect import bisect_left, bisect_right, insort
from collections import deque
import heapq


# ═══════════════════════════ DSU / UNION-FIND ═════════════════════════════════
# WHEN:  connectivity under merges ("are a and b in the same group", "number of components",
#        Kruskal MST, offline "add edges, ask connectivity"), grouping equal/equivalent things.
# MODIFY: keep extra per-component info in arrays indexed by ROOT and merge it inside union().
#   * Component sum/min/max:  s.val[a] += s.val[b]  inside union.
#   * Parity / "enemy" relations (bipartiteness): store parity-to-parent and flip on find
#     (or use 2n nodes: x and x+n as "x" and "not x").
# TIME:  ~O(alpha(n)) per op.
class DSU:
    def __init__(self, n):
        self.p = list(range(n))
        self.sz = [1] * n
        self.comp = n

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]                       # path halving
            x = p[x]
        return x

    def union(self, a, b):                       # returns True if a merge happened
        a, b = self.find(a), self.find(b)
        if a == b:
            return False
        if self.sz[a] < self.sz[b]:
            a, b = b, a
        self.p[b] = a
        self.sz[a] += self.sz[b]
        self.comp -= 1
        return True

    def same(self, a, b):
        return self.find(a) == self.find(b)

    def size(self, x):
        return self.sz[self.find(x)]


# ═══════════════════════════ FENWICK TREE (BIT) ═══════════════════════════════
# WHEN:  point update + prefix/range SUM (or any invertible op: count, xor).
#        Counting inversions, "how many elements < x so far", order statistics (kth).
# MODIFY: values in a[i] -> add(i, v). For counts of values, compress coordinates first.
#   For prefix MAX (only increases) change += to max and drop range().
# TIME:  O(log n).
class BIT:
    def __init__(self, n):
        self.n = n
        self.t = [0] * (n + 1)

    @classmethod
    def from_list(cls, a):                       # O(n) build
        b = cls(len(a))
        for i, v in enumerate(a):
            b.t[i + 1] += v
            j = (i + 1) + ((i + 1) & -(i + 1))
            if j <= b.n:
                b.t[j] += b.t[i + 1]
        return b

    def add(self, i, v):                         # a[i] += v
        i += 1
        t, n = self.t, self.n
        while i <= n:
            t[i] += v
            i += i & -i

    def prefix(self, i):                         # sum of a[0:i]
        s, t = 0, self.t
        while i > 0:
            s += t[i]
            i -= i & -i
        return s

    def range(self, l, r):                       # sum of a[l:r]
        return self.prefix(r) - self.prefix(l)

    def kth(self, k):
        """Smallest 0-based index i with prefix(i+1) >= k (k >= 1). Requires all values >= 0.
        Use with counts: k-th smallest element in a multiset of compressed values."""
        pos, pw = 0, 1 << self.n.bit_length()
        t = self.t
        while pw:
            nxt = pos + pw
            if nxt <= self.n and t[nxt] < k:
                pos = nxt
                k -= t[nxt]
            pw >>= 1
        return pos


class BIT_RangeUpdate:
    """Range add on [l, r), range sum query. WHEN: 'add x to a segment', 'sum of a segment'."""
    def __init__(self, n):
        self.b1, self.b2 = BIT(n + 1), BIT(n + 1)

    def range_add(self, l, r, v):
        self.b1.add(l, v)
        self.b1.add(r, -v)
        self.b2.add(l, v * l)
        self.b2.add(r, -v * r)

    def _pref(self, i):
        return self.b1.prefix(i) * i - self.b2.prefix(i)

    def range_sum(self, l, r):
        return self._pref(r) - self._pref(l)


# ════════════════════ SEGMENT TREE (iterative, any monoid) ════════════════════
# WHEN:  point update + range query for min / max / sum / gcd / xor / "(count, best)" pairs...
# MODIFY: choose op (associative) and identity e.
#     sum:  SegTree(a, lambda x, y: x + y, 0)
#     min:  SegTree(a, min, float('inf'))        max: SegTree(a, max, float('-inf'))
#     gcd:  SegTree(a, math.gcd, 0)
#   For non-trivial data store tuples, e.g. (best_prefix, best_suffix, total, best) for max subarray.
# TIME:  O(log n) per op, O(n) build.
class SegTree:
    def __init__(self, data, op, e):
        self.n = len(data)
        self.op, self.e = op, e
        self.size = 1
        while self.size < self.n:
            self.size *= 2
        self.t = [e] * (2 * self.size)
        self.t[self.size:self.size + self.n] = data
        for i in range(self.size - 1, 0, -1):
            self.t[i] = op(self.t[2 * i], self.t[2 * i + 1])

    def set(self, i, v):
        i += self.size
        self.t[i] = v
        i >>= 1
        while i:
            self.t[i] = self.op(self.t[2 * i], self.t[2 * i + 1])
            i >>= 1

    def get(self, i):
        return self.t[i + self.size]

    def query(self, l, r):                       # op over a[l:r]
        l += self.size
        r += self.size
        a = b = self.e
        while l < r:
            if l & 1:
                a = self.op(a, self.t[l])
                l += 1
            if r & 1:
                r -= 1
                b = self.op(self.t[r], b)
            l >>= 1
            r >>= 1
        return self.op(a, b)


# ═══════════════════════ LAZY SEGMENT TREE (AtCoder-style) ════════════════════
# WHEN:  RANGE UPDATE + RANGE QUERY (range add / range assign / range flip with sum/min/max).
# MODIFY: supply 5 things:
#     op(x, y)           combine two node values        e = identity of op
#     mapping(f, x)      apply lazy tag f to node value x
#     composition(f, g)  tag f applied AFTER tag g  -> single tag
#     id_                identity tag
#   Two ready recipes are given below the class.
class LazySegTree:
    def __init__(self, op, e, mapping, composition, id_, v):
        self.n = len(v)
        self.log = max(1, (self.n - 1).bit_length())
        self.size = 1 << self.log
        self.op, self.e, self.mapping, self.comp, self.id = op, e, mapping, composition, id_
        self.d = [e] * (2 * self.size)
        self.lz = [id_] * self.size
        for i in range(self.n):
            self.d[self.size + i] = v[i]
        for i in range(self.size - 1, 0, -1):
            self._update(i)

    def _update(self, k):
        self.d[k] = self.op(self.d[2 * k], self.d[2 * k + 1])

    def _all_apply(self, k, f):
        self.d[k] = self.mapping(f, self.d[k])
        if k < self.size:
            self.lz[k] = self.comp(f, self.lz[k])

    def _push(self, k):
        self._all_apply(2 * k, self.lz[k])
        self._all_apply(2 * k + 1, self.lz[k])
        self.lz[k] = self.id

    def set(self, p, x):
        p += self.size
        for i in range(self.log, 0, -1):
            self._push(p >> i)
        self.d[p] = x
        for i in range(1, self.log + 1):
            self._update(p >> i)

    def get(self, p):
        p += self.size
        for i in range(self.log, 0, -1):
            self._push(p >> i)
        return self.d[p]

    def prod(self, l, r):                        # query [l, r)
        if l == r:
            return self.e
        l += self.size
        r += self.size
        for i in range(self.log, 0, -1):
            if ((l >> i) << i) != l:
                self._push(l >> i)
            if ((r >> i) << i) != r:
                self._push((r - 1) >> i)
        sl = sr = self.e
        while l < r:
            if l & 1:
                sl = self.op(sl, self.d[l])
                l += 1
            if r & 1:
                r -= 1
                sr = self.op(self.d[r], sr)
            l >>= 1
            r >>= 1
        return self.op(sl, sr)

    def apply(self, l, r, f):                    # apply tag f to [l, r)
        if l == r:
            return
        l += self.size
        r += self.size
        for i in range(self.log, 0, -1):
            if ((l >> i) << i) != l:
                self._push(l >> i)
            if ((r >> i) << i) != r:
                self._push((r - 1) >> i)
        l2, r2 = l, r
        while l < r:
            if l & 1:
                self._all_apply(l, f)
                l += 1
            if r & 1:
                r -= 1
                self._all_apply(r, f)
            l >>= 1
            r >>= 1
        l, r = l2, r2
        for i in range(1, self.log + 1):
            if ((l >> i) << i) != l:
                self._update(l >> i)
            if ((r >> i) << i) != r:
                self._update((r - 1) >> i)


def range_add_range_sum(a):
    """Recipe: node value = (sum, length); tag = add. `t.apply(l,r,x)` adds x; `t.prod(l,r)[0]` is the sum."""
    return LazySegTree(
        op=lambda x, y: (x[0] + y[0], x[1] + y[1]),
        e=(0, 0),
        mapping=lambda f, x: (x[0] + f * x[1], x[1]),
        composition=lambda f, g: f + g,
        id_=0,
        v=[(x, 1) for x in a])


def range_add_range_min(a):
    """Recipe: node value = min; tag = add. For MAX swap min->max and inf->-inf."""
    INF = float("inf")
    return LazySegTree(
        op=min, e=INF,
        mapping=lambda f, x: x + f if x != INF else x,
        composition=lambda f, g: f + g,
        id_=0, v=list(a))


# For range ASSIGN + range sum: tag = None or value; mapping = lambda f,x: x if f is None else (f*x[1], x[1]);
# composition = lambda f,g: g if f is None else f; id_ = None.


# ═══════════════════════════════ SPARSE TABLE ═════════════════════════════════
# WHEN:  STATIC array, many range min / max / gcd queries (idempotent op), no updates.
# TIME:  O(n log n) build, O(1) query.  For range SUM use prefix sums instead.
class SparseTable:
    def __init__(self, a, op=min):
        self.op = op
        self.t = [list(a)]
        j = 1
        while (1 << j) <= len(a):
            prev = self.t[-1]
            half = 1 << (j - 1)
            self.t.append([op(prev[i], prev[i + half]) for i in range(len(a) - (1 << j) + 1)])
            j += 1

    def query(self, l, r):                       # op over a[l:r], l < r
        k = (r - l).bit_length() - 1
        return self.op(self.t[k][l], self.t[k][r - (1 << k)])


# ══════════════════════ MONOTONIC STACK / DEQUE PATTERNS ══════════════════════
# WHEN:  "next greater / smaller element", "span", "largest rectangle in histogram",
#        "sum of subarray minimums", "sliding window max/min".
def next_greater(a):
    """res[i] = index of first j > i with a[j] > a[i], else n."""
    n = len(a)
    res = [n] * n
    st = []
    for i, x in enumerate(a):
        while st and a[st[-1]] < x:
            res[st.pop()] = i
        st.append(i)
    return res


def next_smaller(a):
    """res[i] = index of first j > i with a[j] < a[i], else n."""
    n = len(a)
    res = [n] * n
    st = []
    for i, x in enumerate(a):
        while st and a[st[-1]] > x:
            res[st.pop()] = i
        st.append(i)
    return res


def prev_smaller(a):
    """res[i] = index of last j < i with a[j] < a[i], else -1.  (prev greater: flip the comparison)"""
    n = len(a)
    res = [-1] * n
    st = []
    for i in range(n - 1, -1, -1):
        while st and a[st[-1]] > a[i]:
            res[st.pop()] = i
        st.append(i)
    return res


def largest_rectangle(h):
    """Max area rectangle in a histogram. O(n)."""
    st, best = [], 0
    h = list(h) + [0]
    for i, x in enumerate(h):
        start = i
        while st and st[-1][1] >= x:
            s, hh = st.pop()
            best = max(best, hh * (i - s))
            start = s
        st.append((start, x))
    return best


def sliding_window_max(a, k):
    """max of every window of size k. For MIN flip the comparison. O(n)."""
    dq, res = deque(), []
    for i, x in enumerate(a):
        while dq and a[dq[-1]] <= x:
            dq.pop()
        dq.append(i)
        if dq[0] <= i - k:
            dq.popleft()
        if i >= k - 1:
            res.append(a[dq[0]])
    return res


# ═══════════════════════════════ 2D PREFIX SUMS ═══════════════════════════════
# WHEN:  many rectangle-sum queries on a static grid.
class Prefix2D:
    def __init__(self, g):
        R, C = len(g), len(g[0])
        self.p = [[0] * (C + 1) for _ in range(R + 1)]
        for i in range(R):
            row, run = self.p[i + 1], 0
            for j in range(C):
                run += g[i][j]
                row[j + 1] = self.p[i][j + 1] + run

    def sum(self, r1, c1, r2, c2):               # rows [r1, r2), cols [c1, c2)
        p = self.p
        return p[r2][c2] - p[r1][c2] - p[r2][c1] + p[r1][c1]


# ═══════════════════════════ SORTED LIST (like C++ multiset) ══════════════════
# WHEN:  need insert / delete / "how many elements < x" / "smallest >= x" / k-th, with updates.
#        Python has no built-in balanced BST. (If values are bounded or offline: BIT is faster.)
# TIME:  ~O(sqrt n) per op; fine for n <= ~2e5.
class SortedList:
    def __init__(self, iterable=(), load=500):
        a = sorted(iterable)
        self.load = load
        self.b = [a[i:i + load] for i in range(0, len(a), load)]
        self.mx = [x[-1] for x in self.b]
        self.n = len(a)

    def __len__(self):
        return self.n

    def add(self, x):
        if not self.b:
            self.b, self.mx, self.n = [[x]], [x], 1
            return
        i = bisect_left(self.mx, x)
        if i == len(self.b):
            i -= 1
            self.b[i].append(x)
            self.mx[i] = x
        else:
            insort(self.b[i], x)
        self.n += 1
        if len(self.b[i]) > 2 * self.load:
            half = self.b[i][self.load:]
            del self.b[i][self.load:]
            self.b.insert(i + 1, half)
            self.mx[i] = self.b[i][-1]
            self.mx.insert(i + 1, half[-1])

    def remove(self, x):                         # ValueError if absent
        i = bisect_left(self.mx, x)
        if i == len(self.b):
            raise ValueError(x)
        bk = self.b[i]
        j = bisect_left(bk, x)
        if j == len(bk) or bk[j] != x:
            raise ValueError(x)
        del bk[j]
        self.n -= 1
        if not bk:
            del self.b[i]
            del self.mx[i]
        else:
            self.mx[i] = bk[-1]

    def bisect_left(self, x):                    # number of elements < x
        i = bisect_left(self.mx, x)
        res = sum(len(bk) for bk in self.b[:i])
        return res + (bisect_left(self.b[i], x) if i < len(self.b) else 0)

    def bisect_right(self, x):                   # number of elements <= x
        i = bisect_right(self.mx, x)             # first bucket whose max is > x
        res = sum(len(bk) for bk in self.b[:i])
        return res + (bisect_right(self.b[i], x) if i < len(self.b) else 0)

    def __contains__(self, x):
        i = bisect_left(self.mx, x)
        if i == len(self.b):
            return False
        j = bisect_left(self.b[i], x)
        return j < len(self.b[i]) and self.b[i][j] == x

    def __getitem__(self, k):                    # k-th smallest (0-based); negative ok
        if k < 0:
            k += self.n
        if not 0 <= k < self.n:
            raise IndexError(k)
        for bk in self.b:
            if k < len(bk):
                return bk[k]
            k -= len(bk)

    def pop(self, k=-1):
        x = self[k]
        self.remove(x)
        return x


# ═══════════════════════════ BINARY (XOR) TRIE ════════════════════════════════
# WHEN:  "max XOR of x with any element of a set", XOR of subarrays (insert prefix xors).
# MODIFY: bits = number of bits of the largest value (30 for < 2^30, 31 or 60 as needed).
class XorTrie:
    def __init__(self, bits=30):
        self.B = bits
        self.ch = [[0, 0]]                       # node 0 = root; 0 as a child pointer means "none"
        self.cnt = [0]

    def insert(self, x):
        u = 0
        for b in range(self.B - 1, -1, -1):
            c = (x >> b) & 1
            if not self.ch[u][c]:
                self.ch.append([0, 0])
                self.cnt.append(0)
                self.ch[u][c] = len(self.ch) - 1
            u = self.ch[u][c]
            self.cnt[u] += 1

    def remove(self, x):                         # x must have been inserted
        u = 0
        for b in range(self.B - 1, -1, -1):
            u = self.ch[u][(x >> b) & 1]
            self.cnt[u] -= 1

    def max_xor(self, x):                        # trie must be non-empty
        u, res = 0, 0
        for b in range(self.B - 1, -1, -1):
            c = (x >> b) & 1
            w = self.ch[u][c ^ 1]
            if w and self.cnt[w] > 0:
                res |= 1 << b
                u = w
            else:
                u = self.ch[u][c]
        return res


# ═══════════════════════ HEAPS: lazy deletion, running median ═════════════════
# WHEN:  priority queue where you also need to DELETE arbitrary values (sliding window
#        with heaps, Dijkstra-like with cancellations, "task scheduler" problems).
class RemovableHeap:
    """Min-heap with remove(x). For max-heap push negatives."""
    def __init__(self):
        self.h, self.d, self.size = [], [], 0

    def push(self, x):
        heapq.heappush(self.h, x)
        self.size += 1

    def remove(self, x):                         # x must currently be in the heap
        heapq.heappush(self.d, x)
        self.size -= 1

    def _clean(self):
        while self.d and self.h and self.h[0] == self.d[0]:
            heapq.heappop(self.h)
            heapq.heappop(self.d)

    def top(self):
        self._clean()
        return self.h[0]

    def pop(self):
        self._clean()
        self.size -= 1
        return heapq.heappop(self.h)

    def __len__(self):
        return self.size


class RunningMedian:
    """Insert numbers one by one, ask the median at any time (lower median for even count)."""
    def __init__(self):
        self.lo = []                             # max-heap (negated)
        self.hi = []                             # min-heap

    def add(self, x):
        heapq.heappush(self.lo, -x)
        heapq.heappush(self.hi, -heapq.heappop(self.lo))
        if len(self.hi) > len(self.lo):
            heapq.heappush(self.lo, -heapq.heappop(self.hi))

    def median(self):
        return -self.lo[0]
