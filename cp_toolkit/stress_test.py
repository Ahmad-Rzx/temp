"""
stress_test.py  --  find the smallest failing case for your solution
===============================================================================
HOW TO USE (30 seconds):
  1. Paste your clever solution into fast(case) and a slow-but-obviously-correct one into brute(case).
  2. Write gen() so it returns a RANDOM SMALL case (tiny n, tiny values -> failures are easy to read).
  3. python stress_test.py   -> prints the first case where they disagree.

WHY: when you get Wrong Answer and can't see why, this is faster than staring at code. It also lets
you verify a greedy idea in 2 minutes before you commit to it in the contest.

The demo below checks a (correct) Kadane implementation against an O(n^2) brute force.
"""
import random


def gen():
    n = random.randint(1, 8)
    return [random.randint(-5, 5) for _ in range(n)]


def brute(a):
    best = a[0]
    for i in range(len(a)):
        s = 0
        for j in range(i, len(a)):
            s += a[j]
            best = max(best, s)
    return best


def fast(a):
    best = cur = a[0]
    for x in a[1:]:
        cur = max(x, cur + x)
        best = max(best, cur)
    return best


if __name__ == "__main__":
    random.seed(1)
    for it in range(20000):
        case = gen()
        b, f = brute(list(case)), fast(list(case))      # pass copies if your code mutates the input
        if b != f:
            print("MISMATCH on iteration", it)
            print("input   :", case)
            print("brute   :", b)
            print("fast    :", f)
            break
    else:
        print("all 20000 random cases agree")
