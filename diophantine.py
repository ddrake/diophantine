"""
Module Diophantine

Solve a linear Diophantine equation of the form:

          s*d + t*v = w

where d, v and w are given integers and v > 0
with constraint either p <= s <= q or p <= t <= q.
For a solution to exist, we must have w = m*GCD(d,v) for some integer m.
If this holds, there may be zero or more solutions depending on
the constraint range [p,q] and the coefficient (s or t) that it is applied to.
"""

from math import floor, ceil


def gcd(d, v, qs):
    """
    Compute GCD(d,v) and return it along with the list of quotients.
    where d is an arbitrary integer, v is a positive integer and qs is a list.
    """
    if v <= 0:
        raise ValueError("v must be strictly positive")
    q, r = divmod(d, v)
    if r == 0:
        return v, qs
    qs.append(q)
    return gcd(v, r, qs)


def lincomb(d, v):
    """
    Express GCD(d,v) as a linear combination of d and v
    assuming only that v > 0.  That is, find integers a and b such that
          a*d + b*v = GCD(d,v)
    """
    gd, qs = gcd(d, v, [])
    if len(qs) == 0:
        return gd, 0, 1
    aa = [0, 1]
    bb = [1, -qs[0]]
    for i, q in enumerate(qs[1:]):
        aa.append(aa[i]-q*aa[i+1])
        bb.append(bb[i]-q*bb[i+1])
    return gd, aa[-1], bb[-1]


def diophantine(d, v, w, p, q, v_constraint=False):
    """
    Given a linear Diophantine equation of the form:
             s*d + t*v = w

    along with lower (p) and upper (q) bounds
    on the coefficient of d (or v if v_constraint is True),

    return a list of all (s, t) pairs such that
    a solution exists which satisfies the given constraint.
    """
    if p >= q:
        raise ValueError("p must be less than q")
    sols = []
    gd, _ = gcd(d, v, [])
    m, r = divmod(w, gd)
    if r != 0:
        return []
    gd, a, b = lincomb(d, v)
    if not v_constraint:
        lower = ceil((p-m*a)/v)
        upper = floor((q-m*a)/v)
    elif d > 0:
        lower = ceil((m*b-q)/d)
        upper = floor((m*b-p)/d)
    else:
        lower = ceil((m*b-p)/d)
        upper = floor((m*b-q)/d)

    for n in range(int(lower), int(upper)+1):
        sols.append((m*a+n*v, m*b-n*d))
    return sols


def print_gcd(d, v, gd, a, b):
    print("gcd({0},{1}) = {2} = {3}*{0} + {4}*{1}"
          .format(d, v, gd, a, b))


def print_sol(sol, d, v, w, p, q, v_constraint):
    s, t = sol
    fmt = "(s,t)=({0},{1}) solves {0}*{2} + {1}*{3} = {4}" \
        + " with {5} <= {6} <= {7}"
    print(fmt.format(s, t, d, v, w, p, t if v_constraint else s, q))


if __name__ == '__main__':
    print("GCD Tests")
    dvs = [(237, 141), (6, 3), (5, 2), (3, 7), (199, 98), (98, 199),
           (-129, 273), (-273, 129), (0, 21)]
    for d, v in dvs:
        gd, a, b = lincomb(d, v)
        print_gcd(d, v, gd, a, b)
        assert gd == a*d + b*v

    print()
    print("Solution Tests")
    probs = [(-199, 98, 5, 0, 99, True), (-199, 98, 5, 0, 99, False),
             (23, 31, 4, 0, 99, True), (23, 31, 4, 0, 99, False),
             (23, 31, 4, 0, 3, True), (24, 18, 3, 0, 99, True)]

    for d, v, w, p, q, v_constraint in probs:
        sols = diophantine(d, v, w, p, q, v_constraint)
        if len(sols) == 0:
            print("No solution for s*{} + t*{} = {} with {} <= {} <= {}"
                  .format(d, v, w, p, "t" if v_constraint else "s", q))
        for sol in sols:
            s, t = sol
            msg = "invalid solution"
            assert s*d + t*v == w, msg
            msg = "d coefficient s out of bounds"
            assert p <= s <= q or v_constraint, msg
            msg = "v coefficient t out of bounds"
            assert p <= t <= q or not v_constraint, msg

            print_sol(sol, d, v, w, p, q, v_constraint)
        print()
