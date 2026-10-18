"""
number_theory.py  --  modular arithmetic, primes, factorization, combinatorics
===============================================================================
Every function/class is self-contained: copy only what you need into your solution.
Each block says  WHEN (signals in the problem), MODIFY (what to tweak), TIME.
"""
from math import gcd, isqrt
import random

MOD = 10**9 + 7


# ═════════════════════════════ MODULAR ARITHMETIC ═════════════════════════════
# WHEN:  "print answer modulo 1e9+7 / 998244353", huge powers, divisions under a modulus.
# MODIFY: change MOD. Division a/b mod p  ==  a * mod_inv(b, p) % p  (p prime, b != 0 mod p).
# NOTE:  Python's % is always non-negative for positive modulus, so (a - b) % MOD is safe.
def mod_pow(a, e, m=MOD):
    return pow(a, e, m)                      # built-in is already fast: O(log e)


def mod_inv(a, m=MOD):
    """Inverse of a mod m. Needs gcd(a, m) == 1. (Python 3.8+). For prime m: pow(a, m-2, m)."""
    return pow(a, -1, m)


def inv_table(n, m=MOD):
    """inv[i] = i^-1 mod m for i in 1..n, in O(n). m must be prime and n < m."""
    inv = [0, 1] + [0] * (n - 1)
    for i in range(2, n + 1):
        inv[i] = (m - (m // i) * inv[m % i] % m) % m
    return inv


def ext_gcd(a, b):
    """Returns (g, x, y) with a*x + b*y = g = gcd(a,b).
    WHEN: linear Diophantine equations a*x + b*y = c (solvable iff g | c), inverse for non-prime modulus."""
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def crt(r1, m1, r2, m2):
    """Solve x = r1 (mod m1), x = r2 (mod m2) for ANY moduli.
    Returns (x, lcm) or None if no solution. Chain it for more than two congruences:
        x, m = r[0], mods[0]
        for r_i, m_i in ...: res = crt(x, m, r_i, m_i); if res is None: ...; x, m = res
    """
    g, p, _ = ext_gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    l = m1 // g * m2
    x = (r1 + (r2 - r1) // g * p % (m2 // g) * m1) % l
    return x, l


# ════════════════════════════════ PRIMES ══════════════════════════════════════
# WHEN:  "count / list primes up to n", need primality for many numbers <= ~1e7.
# TIME:  O(n log log n). 1e7 fits in memory using bytearray (10 MB).
def sieve(n):
    """Returns (primes_list, is_prime_bytearray) for 0..n."""
    is_p = bytearray([1]) * (n + 1)
    is_p[0] = 0
    if n >= 1:
        is_p[1] = 0
    for i in range(2, isqrt(n) + 1):
        if is_p[i]:
            is_p[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if is_p[i]], is_p


def spf_sieve(n):
    """spf[x] = smallest prime factor of x. WHEN: factorize MANY numbers <= n quickly (O(log x) each)."""
    spf = list(range(n + 1))
    for i in range(2, isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factorize_spf(x, spf):
    """Prime factorization as {prime: exponent} using an spf table."""
    f = {}
    while x > 1:
        p = spf[x]
        while x % p == 0:
            x //= p
            f[p] = f.get(p, 0) + 1
    return f


def factorize_trial(n):
    """Trial division, O(sqrt n). WHEN: a few numbers up to ~1e12. {prime: exponent}"""
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def is_prime(n):
    """Deterministic Miller-Rabin, correct for all n < 3.3e24.  WHEN: n up to 1e18."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_rho(n):
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    while True:
        c = random.randrange(1, n)
        x = y = random.randrange(2, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d


def factorize_big(n, f=None):
    """Pollard-rho factorization. WHEN: n up to 1e18 and trial division is too slow.
    Returns {prime: exponent}."""
    if f is None:
        f = {}
    if n == 1:
        return f
    if is_prime(n):
        f[n] = f.get(n, 0) + 1
        return f
    d = _pollard_rho(n)
    factorize_big(d, f)
    factorize_big(n // d, f)
    return f


def divisors(n):
    """All divisors of n (sorted), O(sqrt n). For many queries use factorize + product of exponents."""
    small, large = [], []
    i = 1
    while i * i <= n:
        if n % i == 0:
            small.append(i)
            if i != n // i:
                large.append(n // i)
        i += 1
    return small + large[::-1]


def phi_sieve(n):
    """phi[i] = Euler totient for all i <= n.  WHEN: counting coprimes, a^b mod m with huge b."""
    phi = list(range(n + 1))
    for i in range(2, n + 1):
        if phi[i] == i:                          # i is prime
            for j in range(i, n + 1, i):
                phi[j] -= phi[j] // i
    return phi


def euler_phi(n):
    """phi(n) for a single n, O(sqrt n)."""
    res, p = n, 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            res -= res // p
        p += 1
    if n > 1:
        res -= res // n
    return res


def divisor_count_table(n):
    """cnt[i] = number of divisors of i for i <= n. O(n log n). (sum-of-divisors: add j instead of 1)"""
    cnt = [0] * (n + 1)
    for i in range(1, n + 1):
        for j in range(i, n + 1, i):
            cnt[j] += 1
    return cnt


def floor_blocks(n):
    """Yields (l, r, q): for all i in [l, r], n // i == q.  O(sqrt n) blocks.
    WHEN: sum over i of (n // i) * f(i), counting pairs with floor division."""
    l = 1
    while l <= n:
        q = n // l
        r = n // q
        yield l, r, q
        l = r + 1


# ═════════════════════════════ COMBINATORICS ══════════════════════════════════
# WHEN:  "number of ways to choose / arrange ... mod p", n up to ~1e6 and many queries.
# MODIFY: size N must be >= largest n you will ever use (for Catalan: 2*n).
# TIME:  O(N) build, O(1) per query.   mod must be PRIME and > N.
class Comb:
    def __init__(self, N, mod=MOD):
        self.mod = mod
        f = [1] * (N + 1)
        for i in range(1, N + 1):
            f[i] = f[i - 1] * i % mod
        inv = [1] * (N + 1)
        inv[N] = pow(f[N], mod - 2, mod)
        for i in range(N, 0, -1):
            inv[i - 1] = inv[i] * i % mod
        self.f, self.inv = f, inv

    def C(self, n, r):                           # n choose r
        if r < 0 or r > n or n < 0:
            return 0
        return self.f[n] * self.inv[r] % self.mod * self.inv[n - r] % self.mod

    def P(self, n, r):                           # permutations nPr
        if r < 0 or r > n:
            return 0
        return self.f[n] * self.inv[n - r] % self.mod

    def H(self, n, r):                           # multisets: n kinds, choose r with repetition
        return self.C(n + r - 1, r) if n > 0 else int(r == 0)

    def catalan(self, n):                        # needs N >= 2n+1
        return self.C(2 * n, n) * self.f[n] % self.mod * self.inv[n + 1] % self.mod


def lucas(n, r, p, fact):
    """C(n, r) mod SMALL prime p (p ~ <= 1e6) for huge n, r. fact = factorials mod p up to p-1."""
    res = 1
    while n or r:
        ni, ri = n % p, r % p
        if ri > ni:
            return 0
        res = res * fact[ni] % p * pow(fact[ri] * fact[ni - ri] % p, p - 2, p) % p
        n //= p
        r //= p
    return res


def pascal(n, mod=None):
    """Full Pascal triangle C[i][j] for i,j <= n. WHEN: n <= ~2000 and mod is NOT prime."""
    C = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        C[i][0] = 1
        for j in range(1, i + 1):
            v = C[i - 1][j - 1] + C[i - 1][j]
            C[i][j] = v % mod if mod else v
    return C


# ═══════════════════════════ MATRIX EXPONENTIATION ════════════════════════════
# WHEN:  linear recurrence with HUGE n (n up to 1e18): Fibonacci-like, "number of walks of length k".
# MODIFY: build the transition matrix M so that  state_{k+1} = M * state_k, answer = M^n * state_0.
#   Fibonacci: M = [[1,1],[1,0]].  a_n = c1*a_{n-1} + c2*a_{n-2} + c3*a_{n-3}:
#   M = [[c1,c2,c3],[1,0,0],[0,1,0]].
# TIME:  O(k^3 log n) for a k x k matrix.
def mat_mul(A, B, mod=MOD):
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) % mod for col in Bt] for row in A]


def mat_pow(A, e, mod=MOD):
    n = len(A)
    R = [[int(i == j) for j in range(n)] for i in range(n)]
    while e:
        if e & 1:
            R = mat_mul(R, A, mod)
        A = mat_mul(A, A, mod)
        e >>= 1
    return R


def fib(n, mod=MOD):
    """F(n) mod `mod` by fast doubling. O(log n)."""
    def go(k):
        if k == 0:
            return 0, 1
        a, b = go(k >> 1)
        c = a * ((2 * b - a) % mod) % mod
        d = (a * a + b * b) % mod
        return (d, (c + d) % mod) if k & 1 else (c, d)
    return go(n)[0]


# ═══════════════════════════════ MISC HELPERS ═════════════════════════════════
def to_base(n, b):
    """Digits of n in base b (most significant first)."""
    if n == 0:
        return [0]
    out = []
    while n:
        out.append(n % b)
        n //= b
    return out[::-1]


def is_perfect_square(n):
    return n >= 0 and isqrt(n) ** 2 == n


def lcm(a, b):
    return a // gcd(a, b) * b

# QUICK FACTS
#  * gcd(a,b) == gcd(b, a%b).  lcm = a*b/gcd.  gcd of a list: functools.reduce(gcd, lst)
#  * Number of divisors of n = prod(e_i + 1). Max divisors for n <= 1e18 is ~103680.
#  * a^b mod m with HUGE b and gcd(a,m)=1:  pow(a, b % phi(m), m)
#  * Sum 1..n = n(n+1)/2 ; sum of squares = n(n+1)(2n+1)/6
#  * Stars & bars: nonneg integer solutions to x1+..+xk = n is C(n+k-1, k-1).
#  * Inclusion-exclusion over subsets: loop mask in range(1, 1<<k), sign = (-1)**popcount(mask).
