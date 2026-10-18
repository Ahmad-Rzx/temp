# CHEATSHEET — "what do I reach for?"

Read this FIRST on contest day. 90% of problems are solved by recognizing the pattern.

---

## 1. Constraints → intended complexity

CPython does roughly 1–3 × 10^7 simple operations per second; PyPy about 10× more.
Budget your solution against the limit, not against the sample.

| n (input size)      | Intended complexity        | Typical technique                                       |
|---------------------|----------------------------|---------------------------------------------------------|
| n ≤ 10              | O(n!)                      | permutations, brute force, backtracking                 |
| n ≤ 20              | O(2^n · n)                 | bitmask DP, subsets                                     |
| n ≤ 40              | O(2^(n/2))                 | meet in the middle                                      |
| n ≤ 100             | O(n^3)                     | Floyd–Warshall, interval DP                             |
| n ≤ 500             | O(n^3) (tight in Python)   | Floyd–Warshall with list tricks, small DP               |
| n ≤ 2 000–5 000     | O(n^2)                     | 2D DP (LCS, knapsack), pairs — PyPy or tight loops      |
| n ≤ 10^5 – 2·10^5   | O(n log n)                 | sort, heap, BIT/segtree, binary search, Dijkstra        |
| n ≤ 10^6            | O(n) or light O(n log n)   | sieve, prefix sums, two pointers, bulk input            |
| n ≤ 10^9            | O(√n) or O(log n)          | trial division, binary search, fast power               |
| n ≤ 10^18           | O(log n)                   | math formula, matrix power, binary search on answer     |

---

## 2. Problem signal → technique → file

| If the statement says / looks like…                                   | Reach for…                                            | Where |
|-----------------------------------------------------------------------|-------------------------------------------------------|-------|
| "minimum X such that it is possible" / "maximize the minimum"         | **Binary search on the answer** + a `check(x)`        | `search_and_misc.py` |
| sorted order helps, pick best pairs/intervals, "schedule", "assign"   | **Greedy** (sort by the right key), prove w/ exchange | `search_and_misc.py`, `stress_test.py` |
| subarray / substring with sum, distinct count, "longest/shortest"     | **Two pointers / sliding window**                     | `search_and_misc.py` |
| many range-sum queries, no updates                                    | **Prefix sums** (2D: `Prefix2D`)                      | `data_structures.py` |
| range add, one final readout                                          | **Difference array**                                  | `search_and_misc.py` |
| range min/max, static array                                           | **Sparse table**                                      | `data_structures.py` |
| point update + range sum/count                                        | **Fenwick (BIT)**                                     | `data_structures.py` |
| point update + range min/max/gcd, or range update + range query       | **Segment tree / Lazy segment tree**                  | `data_structures.py` |
| "next greater element", histogram, window max                         | **Monotonic stack / deque**                           | `data_structures.py` |
| groups merge, "are they connected", components count                  | **DSU**                                               | `data_structures.py` |
| count inversions, "min adjacent swaps"                                | **Merge sort / BIT**                                  | `search_and_misc.py` |
| "k-th smallest" with inserts/deletes                                  | **BIT.kth / SortedList / two heaps**                  | `data_structures.py` |
| fewest steps, grid, unweighted graph                                  | **BFS** (multi-source for spreading)                  | `graphs.py` |
| weights 0/1                                                           | **0-1 BFS**                                           | `graphs.py` |
| shortest path, non-negative weights                                   | **Dijkstra** (state-space if extra state)             | `graphs.py` |
| negative weights / detect negative cycle                              | **Bellman–Ford**                                      | `graphs.py` |
| all pairs, n ≤ 400                                                    | **Floyd–Warshall**                                    | `graphs.py` |
| "A before B" dependencies, cycle detection in directed graph          | **Topological sort**                                  | `graphs.py` |
| mutually reachable groups                                             | **SCC** (then DP on DAG)                              | `graphs.py` |
| each variable has 2 choices + pairwise conflicts                      | **2-SAT**                                             | `graphs.py` |
| critical edge/vertex whose removal disconnects                        | **Bridges / articulation points**                     | `graphs.py` |
| connect everything cheaply / minimize the max edge                    | **MST (Kruskal)**                                     | `graphs.py` |
| tree: subtree, depth, distances, ancestors                            | **BFS order + reversed-order DP, Euler tour, LCA**    | `graphs.py` |
| assignment / matching / "min cut" / edge-disjoint paths               | **Max flow (Dinic)**, bipartite matching              | `graphs.py` |
| choose items under capacity                                           | **Knapsack** (bitset if only feasibility)             | `dp.py` |
| sequences: increasing / common / transform one into another           | **LIS / LCS / edit distance**                         | `dp.py` |
| n ≤ 20 with "which subset used"                                       | **Bitmask DP**                                        | `dp.py` |
| "count numbers ≤ N with digit property"                               | **Digit DP**                                          | `dp.py` |
| merge/split contiguous segments                                       | **Interval DP**                                       | `dp.py` |
| `dp[i] = min(m*x_i + b)` over earlier j                               | **Li Chao / convex hull trick**                       | `dp.py` |
| answer "mod 1e9+7", counting ways                                     | **Modular math, Comb table, DP with % MOD**           | `number_theory.py` |
| huge n recurrence (1e18)                                              | **Matrix exponentiation**                             | `number_theory.py` |
| primes, divisors, gcd/lcm, coprime counting                           | **Sieve, spf, phi, ext_gcd**                          | `number_theory.py` |
| a few congruences / "x mod a = ..., x mod b = ..."                    | **CRT**                                               | `number_theory.py` |
| find / count a pattern, periods, borders                              | **KMP / Z**                                           | `strings_algo.py` |
| compare substrings fast, "repeated substring of length L"             | **Rolling hash** (+ binary search on L)               | `strings_algo.py` |
| palindromes                                                           | **Manacher**                                          | `strings_algo.py` |
| distinct substrings, longest repeated substring                       | **Suffix array + LCP**                                | `strings_algo.py` |
| many patterns at once                                                 | **Aho–Corasick**                                      | `strings_algo.py` |
| max XOR pair / subarray XOR                                           | **XOR trie + prefix xor**                             | `data_structures.py` |
| points, polygons, hull                                                | **cross product, monotone chain**                     | `geometry.py` |
| many offline range queries with cheap add/remove                      | **Mo's algorithm**                                    | `search_and_misc.py` |
| "n ≤ 40", pick subset with sum target                                 | **Meet in the middle**                                | `search_and_misc.py` |

---

## 3. The 5-minute method for any problem

1. **Restate** in one sentence: input → what exactly to output.
2. **Do the sample by hand**, then invent a tiny case of your own (n = 1, n = 2).
3. **Write the brute force in your head** (what would you do for n ≤ 10?). That is your fallback for partial score.
4. **Look at the constraints** (table 1) → which complexity is required → which technique (table 2).
5. **Edge cases before coding:** n = 1, all equal, already sorted, max values, empty answer (-1?), duplicates.
6. After coding, **run the samples**, then stress-test against brute if anything feels shaky.

---

## 4. Python pitfalls that cost people the contest

| Pitfall | Fix |
|---------|-----|
| `input()` in a loop of 1e5+ lines → TLE | `input = sys.stdin.readline` or read everything at once |
| `print()` inside a loop for 1e5 lines | collect in a list, one `sys.stdout.write("\n".join(...))` |
| `grid = [[0]*m]*n` → all rows are the SAME list | `[[0]*m for _ in range(n)]` |
| recursion depth > ~1000 → RuntimeError | iterative DFS/BFS, or `setrecursionlimit` + `threading.stack_size` |
| `list.pop(0)` / `list.insert(0, x)` is O(n) | `collections.deque` |
| `x in some_list` is O(n) | use `set` / `dict` |
| `s += ch` in a loop is slow-ish | build a list, `"".join` |
| `/` returns float | use `//` for integer division (floats lose precision above 2^53) |
| `-7 // 2 == -4`, `-7 % 3 == 2` (Python floors) | remember it; `int(-7 / 2)` truncates but uses floats |
| Dijkstra without the `if d > dist[u]: continue` check | can blow up to O(E²) |
| Sorting tuples sorts by 1st, then 2nd... | use `key=` for custom orders; `reverse=True` or negate |
| `a.sort()` returns None | `b = sorted(a)` if you need a copy |
| Mutating a list while iterating it | iterate over a copy |
| Default `dict[k]` KeyError | `defaultdict`, `.get(k, 0)`, `Counter` |
| Comparing floats with `==` | use integers, or `abs(a-b) < 1e-9` |
| Forgetting `% MOD` after subtraction | Python gives non-negative mod already: `(a - b) % MOD` is fine |
| Memory: 1e7-element list of ints ≈ 80 MB+ | `bytearray`, `array('i')`, or `numpy` if allowed |
| `lru_cache` recursion on deep DP | convert to bottom-up loop |

**Speed tricks:** local variable aliases inside hot loops (`append = lst.append`), list comprehensions over
for-loops, `sum(map(...))`, `itertools.accumulate`, slice assignment for sieves, Python big ints as bitsets,
`sorted()` + `bisect` instead of hand-written trees, `heapq` for any priority queue.

---

## 5. Contest-day strategy (online round, HackerRank-style scoring)

*Check the actual contest rules page — these are the usual conventions, not guarantees.*

* **Skim ALL problems in the first 5 minutes.** Note which look easiest. Order = easiest first.
* Scoring is usually points per problem (often partial credit per test case) with ties broken by **time**.
  So: fast, correct solves of easy/medium problems beat one heroic hard problem.
* **Partial credit is real.** If you can't find the full solution in ~20 minutes, submit the brute force
  (it often passes the small test cases) and move on. Come back at the end.
* **Timebox:** ~25–30 min per problem before switching. Stuck = write the brute force, submit, move on.
* **Use "Run Code" with the sample** before every submit; compare output format character by character.
* The **time limit shown may differ per language** — read it on the problem page. If a PyPy option exists, test it.
* If a problem is "function completion" style, don't touch their I/O code.
* Don't leave `print` debugging in. Don't print extra spaces/lines unless the format demands it.
* Hard problem last 20 minutes: re-read easy ones for edge cases (n = 1, overflow of ranges, off-by-one).
* Sleep before the contest. A rested brain finds the greedy insight; a tired one rewrites Dijkstra badly.

---

## 6. NOT included (and what to do instead)

Too big to be worth it with two days of prep, and rare in a first-round online contest:
HLD / centroid decomposition, suffix automaton, FFT/NTT (use Python big-int multiplication tricks if
desperate), min-cost max-flow, persistent / wavelet trees, simplex. If a problem needs one of these,
look for the simpler intended trick first — online rounds usually avoid them.
