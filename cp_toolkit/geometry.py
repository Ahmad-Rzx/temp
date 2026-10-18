"""
geometry.py  --  integer-coordinate computational geometry (no floating point errors)
===============================================================================
Points are (x, y) tuples of ints. Prefer integer math (cross/dot products); avoid sqrt and float
compares: compare squared distances instead.
"""
from math import atan2, hypot


def cross(o, a, b):
    """Cross product of OA x OB. > 0: left turn (counter-clockwise), < 0: right turn, 0: collinear."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def dist2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2


def convex_hull(points):
    """Andrew's monotone chain. Returns hull vertices in counter-clockwise order, collinear points removed.
    MODIFY: change `<= 0` to `< 0` to KEEP collinear points on the boundary.
    WHEN: perimeter/area of the hull, 'is point inside all points', farthest pair (rotating calipers on hull)."""
    pts = sorted(set(points))
    if len(pts) <= 2:
        return pts
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def polygon_area2(poly):
    """TWICE the signed area (shoelace) -> an integer. Positive if counter-clockwise."""
    n = len(poly)
    return sum(poly[i][0] * poly[(i + 1) % n][1] - poly[(i + 1) % n][0] * poly[i][1] for i in range(n))


def on_segment(a, b, p):
    return cross(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) \
        and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def segments_intersect(a, b, c, d):
    """Do closed segments ab and cd share at least one point?"""
    d1, d2 = cross(a, b, c), cross(a, b, d)
    d3, d4 = cross(c, d, a), cross(c, d, b)
    if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
        return True
    return on_segment(a, b, c) or on_segment(a, b, d) or on_segment(c, d, a) or on_segment(c, d, b)


def point_in_polygon(poly, p):
    """Returns 1 inside, 0 on the boundary, -1 outside. Works for any simple polygon (ray casting)."""
    n = len(poly)
    inside = False
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        if on_segment(a, b, p):
            return 0
        if (a[1] > p[1]) != (b[1] > p[1]):
            # x-coordinate of the edge at height p[1] compared with p[0], done with integers
            lhs = (p[0] - a[0]) * (b[1] - a[1])
            rhs = (b[0] - a[0]) * (p[1] - a[1])
            if (b[1] > a[1] and lhs < rhs) or (b[1] < a[1] and lhs > rhs):
                inside = not inside
    return 1 if inside else -1


def sort_by_angle(points, center=(0, 0)):
    """Sort points counter-clockwise around center (float atan2; fine unless coordinates are huge)."""
    return sorted(points, key=lambda p: atan2(p[1] - center[1], p[0] - center[0]))


def polygon_perimeter(poly):
    n = len(poly)
    return sum(hypot(poly[i][0] - poly[(i + 1) % n][0], poly[i][1] - poly[(i + 1) % n][1]) for i in range(n))

# FACTS: Pick's theorem: Area = I + B/2 - 1.  Lattice points on segment = gcd(|dx|, |dy|).
#        Closest pair of points: sort by x + divide & conquer, or grid hashing (not included).
