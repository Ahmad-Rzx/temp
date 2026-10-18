"""
strings_algo.py  --  KMP/Z, hashing, Manacher, suffix array + LCP, trie, Aho-Corasick
===============================================================================
"""
import random


# ═════════════════════ PREFIX FUNCTION (KMP) & Z-FUNCTION ═════════════════════
# WHEN:  pattern matching, borders, "smallest period of a string", "count occurrences of each prefix",
#        string compression problems.
# MODIFY: matching text for pattern P -> run on P + '#' + text and look for pi == len(P) (or z == len(P)).
#         smallest period = n - pi[-1] (if it divides n, the string is a repetition).
def prefix_function(s):
    """pi[i] = length of the longest proper prefix of s[:i+1] that is also a suffix of it."""
    n = len(s)
    pi = [0] * n
    for i in range(1, n):
        k = pi[i - 1]
        while k and s[i] != s[k]:
            k = pi[k - 1]
        if s[i] == s[k]:
            k += 1
        pi[i] = k
    return pi


def kmp_find_all(text, pat):
    """Start indices of all occurrences of pat in text. O(n+m)."""
    if not pat:
        return list(range(len(text) + 1))
    pi = prefix_function(pat + "\x00" + text)
    m = len(pat)
    return [i - 2 * m for i in range(2 * m, len(pi)) if pi[i] == m]


def z_function(s):
    """z[i] = length of the longest common prefix of s and s[i:]. z[0] = n."""
    n = len(s)
    z = [0] * n
    l = r = 0
    for i in range(1, n):
        if i < r:
            z[i] = min(r - i, z[i - l])
        while i + z[i] < n and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] > r:
            l, r = i, i + z[i]
    if n:
        z[0] = n
    return z


# ═══════════════════════════ ROLLING HASH ═════════════════════════════════════
# WHEN:  compare substrings in O(1), "is there a repeated substring of length L" (binary search on L),
#        palindromes via forward+reverse hash, counting distinct substrings of one length.
# MODIFY: modulus 2^61-1 with a RANDOM base is safe against anti-hash tests.
class RollingHash:
    MOD = (1 << 61) - 1

    def __init__(self, s, base=None):
        self.B = base or random.randrange(10**6, self.MOD - 1)
        M = self.MOD
        n = len(s)
        self.h = [0] * (n + 1)
        self.pw = [1] * (n + 1)
        for i, ch in enumerate(s):
            c = ch if isinstance(ch, int) else ord(ch)
            self.h[i + 1] = (self.h[i] * self.B + c) % M
            self.pw[i + 1] = self.pw[i] * self.B % M

    def get(self, l, r):                         # hash of s[l:r]
        return (self.h[r] - self.h[l] * self.pw[r - l]) % self.MOD


# ═══════════════════════════════ MANACHER ═════════════════════════════════════
# WHEN:  longest palindromic substring, count ALL palindromic substrings, palindrome partition prep.
def manacher(s):
    """Returns p over the transformed string t = '#a#b#...#'. p[i] = palindrome radius at i,
    which equals the palindrome LENGTH in the original string centered there."""
    t = "#" + "#".join(s) + "#"
    n = len(t)
    p = [0] * n
    c = r = 0
    for i in range(n):
        if i < r:
            p[i] = min(r - i, p[2 * c - i])
        while i - p[i] - 1 >= 0 and i + p[i] + 1 < n and t[i - p[i] - 1] == t[i + p[i] + 1]:
            p[i] += 1
        if i + p[i] > r:
            c, r = i, i + p[i]
    return p


def longest_palindrome(s):
    """Returns the longest palindromic substring."""
    if not s:
        return ""
    p = manacher(s)
    i = max(range(len(p)), key=lambda k: p[k])
    start = (i - p[i]) // 2
    return s[start:start + p[i]]


def count_palindromic_substrings(s):
    """Number of (position, length) palindromic substrings = sum over centers of ceil(p/2)."""
    return sum((x + 1) // 2 for x in manacher(s))


# ═══════════════════════ SUFFIX ARRAY + LCP ═══════════════════════════════════
# WHEN:  count DISTINCT substrings, longest repeated substring, k-th substring, compare suffixes,
#        longest common substring of two strings (join with a separator).
# TIME:  O(n log^2 n) with Python sort: fine to n ~ 1e5 in PyPy, ~5e4 in CPython.
def suffix_array(s):
    n = len(s)
    if n == 0:
        return []
    rank = [ch if isinstance(ch, int) else ord(ch) for ch in s]
    M = max(max(rank), n) + 2
    sa = list(range(n))
    k = 1
    while True:
        key = [rank[i] * M + (rank[i + k] + 1 if i + k < n else 0) for i in range(n)]
        sa.sort(key=key.__getitem__)
        new = [0] * n
        for i in range(1, n):
            new[sa[i]] = new[sa[i - 1]] + (key[sa[i]] != key[sa[i - 1]])
        rank = new
        if rank[sa[-1]] == n - 1:
            break
        k <<= 1
    return sa


def lcp_array(s, sa):
    """Kasai. lcp[i] = LCP(suffix sa[i], suffix sa[i+1]); length n-1."""
    n = len(s)
    rank = [0] * n
    for i, p in enumerate(sa):
        rank[p] = i
    lcp = [0] * max(n - 1, 0)
    k = 0
    for i in range(n):
        if rank[i] == n - 1:
            k = 0
            continue
        j = sa[rank[i] + 1]
        while i + k < n and j + k < n and s[i + k] == s[j + k]:
            k += 1
        lcp[rank[i]] = k
        if k:
            k -= 1
    return lcp


def count_distinct_substrings(s):
    sa = suffix_array(s)
    n = len(s)
    return n * (n + 1) // 2 - sum(lcp_array(s, sa))


# ═══════════════════════════════ TRIE ═════════════════════════════════════════
# WHEN:  many words: prefix counts, "does any word start with ...", longest common prefix, autocomplete.
class Trie:
    def __init__(self):
        self.ch = [{}]
        self.pref = [0]                          # how many inserted words pass through the node
        self.end = [0]                           # how many words end here

    def insert(self, w):
        u = 0
        for c in w:
            nxt = self.ch[u].get(c)
            if nxt is None:
                nxt = len(self.ch)
                self.ch.append({})
                self.pref.append(0)
                self.end.append(0)
                self.ch[u][c] = nxt
            u = nxt
            self.pref[u] += 1
        self.end[u] += 1

    def _walk(self, w):
        u = 0
        for c in w:
            u = self.ch[u].get(c)
            if u is None:
                return None
        return u

    def count_word(self, w):
        u = self._walk(w)
        return self.end[u] if u is not None else 0

    def count_prefix(self, p):
        u = self._walk(p)
        return self.pref[u] if u is not None else 0


# ═══════════════════════════ AHO-CORASICK ═════════════════════════════════════
# WHEN:  find / count occurrences of MANY patterns in a text at once.
# MODIFY: alphabet size A and base char. Total pattern length sets memory (nodes * A).
class AhoCorasick:
    def __init__(self, patterns, alphabet=26, base="a"):
        self.A, self.base = alphabet, ord(base)
        self.nxt = [[-1] * alphabet]
        self.fail = [0]
        self.end_node = []
        for p in patterns:
            u = 0
            for ch in p:
                c = ord(ch) - self.base
                if self.nxt[u][c] == -1:
                    self.nxt[u][c] = len(self.nxt)
                    self.nxt.append([-1] * alphabet)
                    self.fail.append(0)
                u = self.nxt[u][c]
            self.end_node.append(u)
        self.bfs_order = []
        q = []
        for c in range(alphabet):
            v = self.nxt[0][c]
            if v == -1:
                self.nxt[0][c] = 0
            else:
                self.fail[v] = 0
                q.append(v)
        head = 0
        while head < len(q):
            u = q[head]
            head += 1
            self.bfs_order.append(u)
            for c in range(alphabet):
                v = self.nxt[u][c]
                if v == -1:
                    self.nxt[u][c] = self.nxt[self.fail[u]][c]
                else:
                    self.fail[v] = self.nxt[self.fail[u]][c]
                    q.append(v)

    def count_occurrences(self, text):
        """occ[i] = number of occurrences of pattern i in text (overlaps counted)."""
        hits = [0] * len(self.nxt)
        u = 0
        for ch in text:
            u = self.nxt[u][ord(ch) - self.base]
            hits[u] += 1
        for u in reversed(self.bfs_order):
            hits[self.fail[u]] += hits[u]
        return [hits[e] for e in self.end_node]


# ═══════════════════════════ MISC STRING TRICKS ═══════════════════════════════
def min_rotation(s):
    """Index of the lexicographically smallest rotation of s. O(n)."""
    n = len(s)
    s2 = s + s
    i, j, k = 0, 1, 0
    while i < n and j < n and k < n:
        a, b = s2[i + k], s2[j + k]
        if a == b:
            k += 1
            continue
        if a > b:
            i += k + 1
        else:
            j += k + 1
        if i == j:
            j += 1
        k = 0
    return min(i, j)

# QUICK FACTS
#  * Anagram check / grouping: sorted(s) or tuple of 26 counts as a dict key.
#  * Never build strings with += in a loop; collect pieces in a list and "".join(...).
#  * Palindrome check: s == s[::-1].   Counter(s) == Counter(t) for anagram.
#  * ord('a') = 97. Index of letter: ord(c) - 97.
#  * Lexicographically smallest/largest after deletions -> monotonic stack (see data_structures.py).
