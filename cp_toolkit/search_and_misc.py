"""
search_and_misc.py  --  binary search patterns, two pointers, compression, inversions, meet-in-the-middle,
                        bit tricks, Mo's algorithm, greedy/interval patterns, permutations
===============================================================================
"""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from heapq import heapify, heappush, heappop
from math import isqrt


# ═══════════════════ BINARY SEARCH ON THE ANSWER (most important pattern) ═════
# WHEN:  "minimum X such that it is possible", "maximum possible minimum", "smallest time/speed/size",
#        and a check(x) is easy + MONOTONIC (once true, stays true as x grows).
# MODIFY: write pred(x) = "is x enough / feasible?"; pick lo/hi as safe bounds of the answer.
#   Minimize:   first_true(lo, hi, pred)     Maximize: last_true(lo, hi, pred)  (pred True for small x)
def first_true(lo, hi, pred):
    """Smallest x in [lo, hi] with pred(x) True (pred: F F F T T T). Returns hi+1 if none."""
    ans = hi + 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if pred(mid):
            ans = mid
            hi = mid - 1
        else:
            lo = mid + 1
    return ans


def last_true(lo, hi, pred):
    """Largest x in [lo, hi] with pred(x) True (pred: T T T F F F). Returns lo-1 if none."""
    ans = lo - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if pred(mid):
            ans = mid
            lo = mid + 1
        else:
            hi = mid - 1
    return ans


def bisect_real(lo, hi, pred, iters=100):
    """Real-valued: smallest x with pred True (monotone F..T). 100 iterations is plenty for doubles."""
    for _ in range(iters):
        mid = (lo + hi) / 2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return hi


def ternary_min_int(lo, hi, f):
    """Minimum of a strictly unimodal (down then up) f on integers [lo, hi]. For max negate f."""
    while hi - lo > 2:
        m1 = lo + (hi - lo) // 3
        m2 = hi - (hi - lo) // 3
        if f(m1) < f(m2):
            hi = m2
        else:
            lo = m1
    return min(range(lo, hi + 1), key=f)


# bisect cheat sheet (sorted list a):
#   bisect_left(a, x)  = number of elements <  x  = first index with a[i] >= x
#   bisect_right(a, x) = number of elements <= x  = first index with a[i] >  x
#   count of x in a = bisect_right(a,x) - bisect_left(a,x)
#   largest element <= x: i = bisect_right(a,x)-1 (check i>=0);  smallest >= x: i = bisect_left(a,x) (check i<len)


# ═════════════════════════ TWO POINTERS / SLIDING WINDOW ══════════════════════
# WHEN:  contiguous subarray/substring with a monotone property (adding extends, removing shrinks):
#        "longest/shortest subarray with sum <= K / >= K", "at most K distinct", "pairs in sorted array".
# MODIFY: keep window state (sum, counter); expand r each step, shrink l while the window is invalid.
def longest_subarray_sum_at_most(a, K):
    """Longest subarray with sum <= K (non-negative elements)."""
    l = s = best = 0
    for r, x in enumerate(a):
        s += x
        while s > K:
            s -= a[l]
            l += 1
        best = max(best, r - l + 1)
    return best


def count_subarrays_at_most_k_distinct(a, k):
    """Number of subarrays with at most k distinct values.  exactly k = f(k) - f(k-1)."""
    cnt = defaultdict(int)
    l = res = 0
    for r, x in enumerate(a):
        cnt[x] += 1
        while len(cnt) > k:
            cnt[a[l]] -= 1
            if cnt[a[l]] == 0:
                del cnt[a[l]]
            l += 1
        res += r - l + 1
    return res


def count_subarrays_sum_equals(a, K):
    """Subarrays with sum == K, works with NEGATIVES (prefix sums + hashmap)."""
    seen = defaultdict(int)
    seen[0] = 1
    s = res = 0
    for x in a:
        s += x
        res += seen[s - K]
        seen[s] += 1
    return res


def pairs_sum_at_most(a, K):
    """Number of index pairs i<j with a[i]+a[j] <= K (order of pairs irrelevant -> sort first)."""
    a = sorted(a)
    l, r, res = 0, len(a) - 1, 0
    while l < r:
        if a[l] + a[r] <= K:
            res += r - l
            l += 1
        else:
            r -= 1
    return res


# ═════════════════════ PREFIX SUMS / DIFFERENCE ARRAY / COMPRESSION ═══════════
# prefix:  pre[i] = a[0]+...+a[i-1];  sum(l..r-1) = pre[r]-pre[l]       (itertools.accumulate)
# difference array: to add v on [l, r) do d[l]+=v; d[r]-=v; then accumulate(d) gives final values.
# xor-prefix works the same: xor of a[l:r] = px[r] ^ px[l].
def compress(a):
    """Coordinate compression. Returns (ranks, sorted_unique): ranks[i] in [0, len(unique)).
    WHEN: values up to 1e9 but you need array/BIT indices."""
    u = sorted(set(a))
    idx = {v: i for i, v in enumerate(u)}
    return [idx[x] for x in a], u


def count_inversions(a):
    """Number of pairs i<j with a[i] > a[j]. O(n log n) merge sort.
    WHEN: 'minimum adjacent swaps to sort' (= inversions), 'how unsorted'. (BIT alternative in data_structures.)"""
    def sort(arr):
        n = len(arr)
        if n <= 1:
            return arr, 0
        m = n // 2
        L, x = sort(arr[:m])
        R, y = sort(arr[m:])
        inv = x + y
        i = j = 0
        out = []
        while i < len(L) and j < len(R):
            if L[i] <= R[j]:
                out.append(L[i])
                i += 1
            else:
                out.append(R[j])
                j += 1
                inv += len(L) - i
        out += L[i:]
        out += R[j:]
        return out, inv
    return sort(list(a))[1]


# ═══════════════════ MEET IN THE MIDDLE / BIT TRICKS ══════════════════════════
# WHEN:  n <= ~40 and 2^n is too much but 2^(n/2) is fine: subset sums, "closest to target".
def subset_sums(arr):
    sums = [0]
    for x in arr:
        sums += [s + x for s in sums]
    return sums


def count_subsets_with_sum(arr, target):
    half = len(arr) // 2
    left, right = subset_sums(arr[:half]), subset_sums(arr[half:])
    cnt = defaultdict(int)
    for s in right:
        cnt[s] += 1
    return sum(cnt[target - s] for s in left)


def submasks(mask):
    """All submasks of mask (including 0 and mask). Total over all masks = 3^n."""
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub - 1) & mask


def popcount(x):
    return bin(x).count("1")                     # int.bit_count() on Python 3.10+

# BIT FACTS:  x & -x = lowest set bit | x & (x-1) clears it | (x >> i) & 1 tests bit i
#             1 << n = 2^n | x ^ y toggles | all-ones of n bits = (1 << n) - 1
#             XOR tricks: x ^ x = 0, find the unique element = xor of all.
#             Subset enumeration: for mask in range(1 << n): for i in range(n): if mask >> i & 1: ...


# ══════════════════════════ MO'S ALGORITHM (offline range queries) ════════════
# WHEN:  many range queries [l, r] where moving an endpoint by one is cheap and no update structure
#        fits: distinct count in range, mode, sum of f(count) ...  Time O((n + q) sqrt n): slow in
#        CPython, ok for n,q <= ~5e4 (PyPy ~2e5).
# MODIFY: provide add(i), remove(i) that update your global state, and get() that returns the answer.
def mo_algorithm(n, queries, add, remove, get):
    """queries: [(l, r)] inclusive 0-indexed. Returns answers in input order."""
    q = len(queries)
    B = max(1, int(n / max(1.0, q ** 0.5)))
    order = sorted(range(q), key=lambda i: (queries[i][0] // B,
                   queries[i][1] if (queries[i][0] // B) % 2 == 0 else -queries[i][1]))
    ans = [None] * q
    cl, cr = 0, -1
    for qi in order:
        l, r = queries[qi]
        while cl > l:
            cl -= 1
            add(cl)
        while cr < r:
            cr += 1
            add(cr)
        while cl < l:
            remove(cl)
            cl += 1
        while cr > r:
            remove(cr)
            cr -= 1
        ans[qi] = get()
    return ans


# ═════════════════════ GREEDY & INTERVAL PATTERNS ═════════════════════════════
# Greedy checklist: sort by the right key (end time, ratio, deadline), then take the best feasible
# choice. If you can't PROVE it (exchange argument), test it against brute force with stress_test.py.
def merge_intervals(iv):
    """Merge overlapping closed intervals [(l, r)]. Touching intervals (r == next l) are merged."""
    out = []
    for l, r in sorted(iv):
        if out and l <= out[-1][1]:
            out[-1][1] = max(out[-1][1], r)
        else:
            out.append([l, r])
    return [tuple(x) for x in out]


def max_non_overlapping(iv):
    """Activity selection: max number of pairwise non-overlapping intervals (half-open [l, r))."""
    iv = sorted(iv, key=lambda x: x[1])
    cnt, last = 0, float("-inf")
    for l, r in iv:
        if l >= last:
            cnt += 1
            last = r
    return cnt


def max_overlap(iv):
    """Maximum number of intervals covering one point = minimum rooms needed. Half-open [l, r).
    Sweep line: +1 at l, -1 at r (process -1 before +1 at equal coordinates)."""
    ev = []
    for l, r in iv:
        ev.append((l, 1))
        ev.append((r, -1))
    ev.sort()
    cur = best = 0
    for _, d in ev:
        cur += d
        best = max(best, cur)
    return best


def min_cost_merge_all(a):
    """Repeatedly merge two items at cost = their sum (Huffman). Min total cost."""
    h = list(a)
    heapify(h)
    cost = 0
    while len(h) > 1:
        x, y = heappop(h), heappop(h)
        cost += x + y
        heappush(h, x + y)
    return cost


def next_permutation(a):
    """In-place next lexicographic permutation of list a. Returns False if a was the last one."""
    i = len(a) - 2
    while i >= 0 and a[i] >= a[i + 1]:
        i -= 1
    if i < 0:
        a.reverse()
        return False
    j = len(a) - 1
    while a[j] <= a[i]:
        j -= 1
    a[i], a[j] = a[j], a[i]
    a[i + 1:] = reversed(a[i + 1:])
    return True

# CUSTOM SORT ORDER:  from functools import cmp_to_key
#   a.sort(key=cmp_to_key(lambda x, y: -1 if x + y > y + x else 1))   # e.g. "largest number" concatenation
# Sort by several keys: a.sort(key=lambda p: (-p[1], p[0]))   (descending score, then ascending name)
