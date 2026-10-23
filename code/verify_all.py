#!/usr/bin/env python3
"""Single verifier for "Nonreal poles of Weniger's delta transformation".

Usage:
    python3 -B -O verify_all.py OUTPUT_DIRECTORY

Writes OUTPUT_DIRECTORY/verification.json (sorted keys, no timestamps, no
absolute paths, floats rounded to 8 significant figures) and prints a short
summary. Every exact check uses Python integers and fractions.Fraction and
raises RuntimeError on failure (no assert statements, so -O is safe).

Sections whose key starts with "numerical_" use mpmath floating-point
arithmetic. They illustrate the practical discussion in the paper and are
not part of any proof. All other sections are exact.
"""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction as F
from functools import reduce
from math import comb, factorial, gcd, lcm
from pathlib import Path


# ---------------------------------------------------------------------------
# Basic exact helpers
# ---------------------------------------------------------------------------

def check(condition, message):
    if not condition:
        raise RuntimeError("check failed: " + message)


def trim(p):
    p = [F(c) for c in p]
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [F(0)]


def add(*polys):
    out = [F(0)] * max(map(len, polys))
    for p in polys:
        for j, c in enumerate(p):
            out[j] += c
    return trim(out)


def scale(p, a):
    return trim([a * c for c in p])


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            out[i + j] += a * b
    return trim(out)


def power(p, n):
    out = [F(1)]
    for _ in range(n):
        out = mul(out, p)
    return out


def evaluate(p, x):
    v = F(0)
    for c in reversed(p):
        v = v * x + c
    return v


def divmod_poly(p, q):
    p, q = trim(p), trim(q)
    check(q != [0], "division by zero polynomial")
    quotient = [F(0)] * max(1, len(p) - len(q) + 1)
    while p != [0] and len(p) >= len(q):
        d = len(p) - len(q)
        a = p[-1] / q[-1]
        quotient[d] = a
        p = add(p, [F(0)] * d + scale(q, -a))
    return trim(quotient), trim(p)


def gcd_poly(p, q):
    while trim(q) != [0]:
        p, q = q, divmod_poly(p, q)[1]
    return scale(p, 1 / trim(p)[-1])


def primitive(p):
    p = trim(p)
    den = lcm(*(c.denominator for c in p))
    ints = [int(c * den) for c in p]
    g = reduce(gcd, (abs(c) for c in ints))
    if ints[-1] < 0:
        g = -g
    return [c // g for c in ints], F(g, den)


def rising(a, k):
    out = F(1)
    for j in range(k):
        out *= a + j
    return out


def determinant(matrix):
    a = [[F(c) for c in row] for row in matrix]
    out = F(1)
    for j in range(len(a)):
        piv = next((i for i in range(j, len(a)) if a[i][j]), None)
        if piv is None:
            return F(0)
        if piv != j:
            a[piv], a[j] = a[j], a[piv]
            out = -out
        p = a[j][j]
        out *= p
        for i in range(j + 1, len(a)):
            m = a[i][j] / p
            for k in range(j + 1, len(a)):
                a[i][k] -= m * a[j][k]
            a[i][j] = 0
    return out


def sylvester(p, q):
    """Sylvester matrix, descending coefficients in each row; rows of p first."""
    m, n = len(p) - 1, len(q) - 1
    p, q = list(reversed(p)), list(reversed(q))
    rows = [[0] * j + p + [0] * (n - 1 - j) for j in range(n)]
    rows += [[0] * j + q + [0] * (m - 1 - j) for j in range(m)]
    return rows


def discriminant_cubic(p):
    d, c, b2, a = p
    return b2 * b2 * c * c - 4 * a * c ** 3 - 4 * b2 ** 3 * d - 27 * a * a * d * d + 18 * a * b2 * c * d


def r8(x):
    """Round a real number to 8 significant figures for deterministic JSON."""
    return float(f"{float(x):.8g}")


# ---------------------------------------------------------------------------
# Exact Sturm counts and gcd over Z[z] (primitive pseudo-remainder sequences)
# ---------------------------------------------------------------------------

def iprim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    g = reduce(gcd, (abs(c) for c in p))
    if g == 0:
        return [0]
    return [c // g for c in p]


def iprem(a, b):
    """Remainder of |lc(b)|^(deg a - deg b + 1) * a modulo b (sign preserving)."""
    a = list(a)
    lb, db = b[-1], len(b) - 1
    delta = len(a) - 1 - db
    if delta < 0:
        return a
    labs, sgn = abs(lb), (1 if lb > 0 else -1)
    for _ in range(delta + 1):
        if len(a) - 1 < db:
            a = [c * labs for c in a]
            continue
        la, shift = a[-1], len(a) - 1 - db
        a = [c * labs for c in a]
        for i, c in enumerate(b):
            a[i + shift] -= sgn * la * c
        a.pop()
        while len(a) > 1 and a[-1] == 0:
            a.pop()
    return a or [0]


def int_poly(p):
    return primitive(p)[0]


def real_root_count(p):
    """(number of distinct real zeros, squarefree flag) by an exact Sturm sequence."""
    p = iprim(int_poly(p))
    seq = [p, iprim([j * p[j] for j in range(1, len(p))])]
    while len(seq[-1]) > 1:
        r = [-c for c in iprem(seq[-2], seq[-1])]
        if all(c == 0 for c in r):
            break
        seq.append(iprim(r))
    at_pos = [q[-1] for q in seq]
    at_neg = [q[-1] * (-1) ** (len(q) - 1) for q in seq]

    def changes(v):
        s = [x for x in v if x != 0]
        return sum(1 for x, y in zip(s, s[1:]) if (x > 0) != (y > 0))

    return changes(at_neg) - changes(at_pos), len(seq[-1]) == 1


def nonreal_count(p):
    """Number of nonreal zeros (with multiplicity) of a squarefree real polynomial."""
    nreal, squarefree = real_root_count(p)
    check(squarefree, "Sturm count requires a squarefree polynomial")
    return len(trim(p)) - 1 - nreal


def distinct_nonreal_count(p):
    """Number of distinct nonreal zeros (works without squarefreeness)."""
    p = trim(p)
    g = gcd_poly(p, [j * p[j] for j in range(1, len(p))])
    sqf, r = divmod_poly(p, g)
    check(r == [0], "squarefree part")
    nreal, ok = real_root_count(sqf)
    check(ok, "squarefree part is squarefree")
    return len(sqf) - 1 - nreal


def gcd_degree(p, q):
    a, b = iprim(int_poly(p)), iprim(int_poly(q))
    if len(a) < len(b):
        a, b = b, a
    while not (len(b) == 1 and b[0] == 0):
        r = iprem(a, b)
        a, b = b, (iprim(r) if any(r) else [0])
    return len(a) - 1


# ---------------------------------------------------------------------------
# The delta transformation
# ---------------------------------------------------------------------------

def delta_polynomials(moments, n, beta, k):
    """P, Q with delta_k^(n)(beta; z) = P(z) / (z^(n+1) Q(z)), from eq. (rational)."""
    q = [F(comb(k, j)) * rising(F(n) + beta + j, k - 1) / moments[n + j + 1] for j in range(k + 1)]
    p = [F(0)] * (n + k + 1)
    for j, a in enumerate(q):
        for ell in range(n + j + 1):
            p[n + j - ell] += a * (-1) ** ell * moments[ell]
    return trim(p), trim(q)


def borghi_P(moments, n, k):
    """Borghi's auxiliary polynomial P_k(z) = sum_j C(k,j) z^j / mu_(n+1+j)."""
    return [F(comb(k, j)) / moments[n + 1 + j] for j in range(k + 1)]


def direct_delta(moments, n, beta, k, z):
    """delta_k^(n)(beta; z) from the finite-difference definition, eq. (delta)."""
    sums = [sum((-1) ** ell * moments[ell] / z ** (ell + 1) for ell in range(n + j + 1)) for j in range(k + 1)]
    diffs = [(-1) ** (n + j + 1) * moments[n + j + 1] / z ** (n + j + 2) for j in range(k + 1)]
    weights = [(-1) ** (k - j) * comb(k, j) * rising(F(n) + beta + j, k - 1) for j in range(k + 1)]
    num = sum(w * s / d for w, s, d in zip(weights, sums, diffs))
    den = sum(w / d for w, d in zip(weights, diffs))
    return num / den


def series_at_infinity(p, q, n, last):
    numerator, denominator = list(reversed(p)), list(reversed(q))
    check(len(p) - 1 == n + len(q) - 1, "degree pattern for series expansion")
    out = []
    for ell in range(last + 1):
        v = numerator[ell] if ell < len(numerator) else F(0)
        for j in range(1, min(ell, len(denominator) - 1) + 1):
            v -= denominator[j] * out[ell - j]
        out.append(v / denominator[0])
    return out


CUBIC_MEASURE = "delta_1 + 50 delta_(1/50) + (1/1000) exp(-t) dt"
QUARTIC_MEASURE = "delta_1 + 100000 delta_(1/100) + (1/1000) exp(-t) dt"


def cubic_moments(N):
    return [F(1) + F(50) ** (1 - ell) + F(factorial(ell), 1000) for ell in range(N)]


def quartic_moments(N):
    return [F(1) + F(10) ** (5 - 2 * ell) + F(factorial(ell), 1000) for ell in range(N)]


# ---------------------------------------------------------------------------
# Section A: symbolic identities used in the proofs
# ---------------------------------------------------------------------------

def mpoly_mul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, F(0)) + ca * cb
    return {e: c for e, c in out.items() if c != 0}


def mpoly_add(*ps):
    out = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, F(0)) + c
    return {e: c for e, c in out.items() if c != 0}


def mpoly_scale(p, s):
    return {e: c * s for e, c in p.items() if c * s != 0}


def mpoly_from_b(poly_b):
    """Univariate polynomial in b -> trivariate dict in (b, x, y)."""
    return {(i, 0, 0): F(c) for i, c in enumerate(poly_b) if c != 0}


def check_cubic_identities():
    D, C, B, A = [0, 1, 1], [6, 9, 3], [18, 15, 3], [12, 7, 1]
    # Corner values of F.
    f00 = mul(power(B, 2), power(C, 2))
    f01 = add(f00, scale(mul(A, power(C, 3)), -4))
    f10 = add(f00, scale(mul(power(B, 3), D), -4))
    f11 = add(f01, scale(mul(power(B, 3), D), -4), scale(mul(mul(A, B), mul(C, D)), 18),
              scale(mul(power(A, 2), power(D, 2)), -27))
    check(f01 == scale(mul(power(C, 2), mul(mul([2, 1], [3, 1]), [2, -5, -1])), 3), "F(0,1)")
    check(f10 == scale(mul(power(B, 2), mul(mul([1, 1], [2, 1]), [6, -3, -1])), 3), "F(1,0)")
    check(f11 == scale(mul([1, 1], [3, 1]), 432), "F(1,1)")
    middle = scale(mul([1, 1], [3, 1]), 2)
    check(add(middle, D) == trim(C) and add(A, middle) == trim(B), "factorisation at x=y=1")
    # NEW: eq. (disc): Disc_zeta(D + C zeta + B x zeta^2 + A x^2 y zeta^3) = x^2 F(x, y) over Q[b, x, y].
    d = mpoly_from_b(D)
    c = mpoly_from_b(C)
    b2 = mpoly_mul(mpoly_from_b(B), {(0, 1, 0): F(1)})
    a = mpoly_mul(mpoly_from_b(A), {(0, 2, 1): F(1)})
    m = mpoly_mul
    disc = mpoly_add(m(m(b2, b2), m(c, c)), mpoly_scale(m(a, m(c, m(c, c))), -4),
                     mpoly_scale(m(b2, m(b2, m(b2, d))), -4), mpoly_scale(m(m(a, a), m(d, d)), -27),
                     mpoly_scale(m(m(a, b2), m(c, d)), 18))
    Bp, Cp, Dp, Ap = (mpoly_from_b(v) for v in (B, C, D, A))
    y = {(0, 0, 1): F(1)}
    x = {(0, 1, 0): F(1)}
    Fxy = mpoly_add(m(m(Bp, Bp), m(Cp, Cp)), mpoly_scale(m(Ap, m(Cp, m(Cp, m(Cp, y)))), -4),
                    m(x, mpoly_add(mpoly_scale(m(Bp, m(Bp, m(Bp, Dp))), -4),
                                   mpoly_scale(m(m(Ap, Bp), m(Cp, m(Dp, y))), 18))),
                    mpoly_scale(m(m(Ap, Ap), m(m(Dp, Dp), m(m(x, x), m(y, y)))), -27))
    check(disc == m(m(x, x), Fxy), "eq. (disc): Disc(p) = x^2 F(x,y) over Q[b,x,y]")
    # Third finite difference of (b+j)_2 vanishes.
    fd = add(*[scale(mul([j, 1], [j + 1, 1]), (-1) ** j * comb(3, j)) for j in range(4)])
    check(fd == [0], "third difference of rising factorial")
    return {
        "corner_identities_F00_F01_F10_F11": "verified over Q[b]",
        "F11_coefficients_ascending": [int(v) for v in f11],
        "disc_identity_Disc_p_equals_x2_F": "verified as a polynomial identity over Q[b,x,y]",
        "disc_identity_number_of_monomials": len(disc),
        "factorisation_at_x_y_one": "verified over Q[b]",
        "positive_boundary": "(sqrt(33)-5)/2, positive root of b^2+5b-2",
    }


def check_deltaQ_and_tailratio():
    """Eq. (DeltaQ) over Q[b] for generic rational moments; eq. (tailratio) over Q[b]."""
    checked = []
    for n in range(0, 4):
        for k in range(1, 9):
            inv = [F(1, 3 + 2 * i * i + i) for i in range(n + k + 3)]  # generic positive 1/mu values
            # LHS coefficient of z^(n+2+j): (-1)^(k-j) C(k,j) (b+j)_(k-1) (-1)^(n+j+1) inv[n+j+1]
            for j in range(k + 1):
                rb = [F(1)]
                for i in range(k - 1):
                    rb = mul(rb, [F(j + i), F(1)])
                lhs = scale(rb, F((-1) ** (k - j) * comb(k, j) * (-1) ** (n + j + 1)) * inv[n + j + 1])
                rhs = scale(rb, F((-1) ** (n + k + 1) * comb(k, j)) * inv[n + j + 1])
                check(lhs == rhs, "eq. (DeltaQ) coefficient")
            checked.append([n, k])
    tails = []
    for k in range(3, 13):
        def h(j):
            out = [F(comb(k, j))]
            for i in range(k - 1):
                out = mul(out, [F(j + i), F(1)])
            return out
        # R_k(b) = 2k/(k-1) * (b+k-1)(b+2k-3)/((b+k-2)(b+2k-2)) = 2k/(k-1) + 2k/((b+k-2)(b+2k-2))
        lhs = mul(mul(power(h(k - 1), 2), [F(k - 2), F(1)]), mul([F(2 * k - 2), F(1)], [F(k - 1)]))
        rhs = mul(mul(h(k - 2), h(k)), scale(mul([F(k - 1), F(1)], [F(2 * k - 3), F(1)]), F(2 * k)))
        check(lhs == rhs, "eq. (tailratio) product form over Q[b]")
        pf = add(scale(mul([F(k - 2), F(1)], [F(2 * k - 2), F(1)]), F(2 * k, k - 1)), [F(2 * k)])
        check(pf == scale(mul([F(k - 1), F(1)], [F(2 * k - 3), F(1)]), F(2 * k, k - 1)),
              "eq. (tailratio) partial-fraction form")
        r0 = F(k * (2 * k - 3), (k - 1) * (k - 2))
        if k >= 4:
            check(r0 <= F(10, 3) < 4, "tail bound")
        tails.append({"k": k, "R_k_at_b0": str(r0)})
    # Middle inequality of eq. (tailbound): 10(k-1)(k-2) - 3k(2k-3) = (k-4)(4k-5) over Q[k].
    check(add(scale(mul([F(-1), F(1)], [F(-2), F(1)]), F(10)), scale(mul([F(0), F(1)], [F(-3), F(2)]), F(-3)))
          == mul([F(-4), F(1)], [F(-5), F(4)]), "tail bound factorisation (k-4)(4k-5)")
    # Rational form (eq. rational) against the raw definition at several points and orders.
    mom = cubic_moments(12)
    pts = 0
    for k in (3, 4, 5):
        for n in (0, 1, 2):
            for beta in (F(1), F(1, 3), F(5, 2)):
                p, q = delta_polynomials(mom, n, beta, k)
                for z in (F(1), F(2), F(7, 3)):
                    check(direct_delta(mom, n, beta, k, z) == evaluate(p, z) / (z ** (n + 1) * evaluate(q, z)),
                          "rational form vs finite-difference definition")
                    pts += 1
    return {"DeltaQ_identity_checked_n_k": checked, "tail_ratio_identity_over_Qb": tails,
            "rational_form_vs_definition_point_checks": pts}


# ---------------------------------------------------------------------------
# Section B: the explicit cubic and quartic examples (ported from the first version)
# ---------------------------------------------------------------------------

def exact_cubic_example():
    moments = cubic_moments(9)
    check(moments[:5] == [F(51001, 1000), F(2001, 1000), F(511, 500), F(629, 625), F(128001, 125000)], "moments")
    p, q = delta_polynomials(moments, 0, F(1), 3)
    pi, ps = primitive(p)
    qi, qs = primitive(q)
    check(qi == [54855937892, 966634623774, 1963236457665, 1071932365000], "Q integer")
    check(pi == [1791106733855748, 46466411172340109, 97982085915007665, 54669622547365000], "P integer")
    check(ps / qs == F(1, 1000), "factor 1000")
    disc = discriminant_cubic(qi)
    check(disc == -16407441394630862189702708459297984636170197900, "cubic discriminant")
    check(gcd_poly(pi, qi) == [1], "gcd")
    res = determinant(sylvester(pi, qi))
    check(res != 0 and res.denominator == 1, "resultant")
    check(series_at_infinity(p, q, 0, 4) == [(-1) ** l * moments[l] for l in range(5)], "accuracy")
    for z in (F(1), F(2), F(7, 3)):
        check(direct_delta(moments, 0, F(1), 3, z) == evaluate(pi, z) / (1000 * z * evaluate(qi, z)), "direct")
    pb = borghi_P(moments, 0, 3)
    check(nonreal_count(pb) == 2, "Borghi P_3 has nonreal zeros (scenario ii)")
    return {"measure": CUBIC_MEASURE, "n": 0, "beta": "1", "k": 3,
            "P_integer_ascending": pi, "Q_integer_ascending": qi,
            "rational_approximant": "P(z)/(1000*z*Q(z))", "discriminant_Q": str(disc),
            "gcd_P_Q": "1", "resultant_P_Q": str(res.numerator),
            "borghi_P3_nonreal_zero_count": 2}


def exact_quartic_example():
    moments = quartic_moments(8)
    p, q = delta_polynomials(moments, 0, F(1, 3), 4)
    pi, ps = primitive(p)
    qi, qs = primitive(q)
    check(qi == [139692321001730, 508388122214334600, 24653951028874577475, 40536707773004833728,
                 17621547662058132000], "Q4 integer")
    check(pi == [-83485463338713495208696, 8862395637097892836188527, 808345456247443334317918249,
                 1345357389237898766149476576, 587390801791674306730044000], "P4 integer")
    check(ps / qs == F(3, 1000), "scale 3/1000 in delta_4 = 3 P/(1000 z Q)")
    deriv = [j * qi[j] for j in range(1, 5)]
    disc = determinant(sylvester(qi, deriv)) / qi[-1]
    check(disc < 0 and disc.denominator == 1, "quartic discriminant negative")
    lag = evaluate(deriv, F(-1)) ** 2 - evaluate(qi, F(-1)) * evaluate([j * deriv[j] for j in range(1, 4)], F(-1))
    check(lag == -16188455819862407838029037564016194858, "Laguerre value")
    check(gcd_poly(pi, qi) == [1], "gcd quartic")
    check(series_at_infinity(p, q, 0, 5) == [(-1) ** l * moments[l] for l in range(6)], "accuracy quartic")
    _, q3 = delta_polynomials(moments, 0, F(1, 3), 3)
    check(discriminant_cubic(q3) > 0, "same measure cubic is real-rooted")
    check(nonreal_count(borghi_P(moments, 0, 4)) == 2, "Borghi P_4 has nonreal zeros (scenario ii)")
    # Elimination polynomial H of the exceptional-moment lemma.
    n, k, N = 0, 4, 4
    H = [F(0)] * k
    for j in range(k):
        for ell in range(n + j + 1, N + 1):
            H[N + j - ell] -= q[j] * (-1) ** ell * moments[ell]
    check(H[0] == -q[0] * (-1) ** N * moments[N] and H[0] != 0 and len(trim(H)) - 1 <= k - 1, "H")
    return {"measure": QUARTIC_MEASURE, "n": 0, "beta": "1/3", "k": 4,
            "P_integer_ascending": pi, "Q_integer_ascending": qi,
            "rational_approximant": f"({ps / qs})*P(z)/(z*Q(z))",
            "discriminant_Q": str(disc.numerator), "laguerre_value_at_minus_one": str(lag),
            "gcd_P_Q": "1", "same_measure_cubic_discriminant_sign": 1,
            "borghi_P4_nonreal_zero_count": 2}


def check_legendre():
    out = []
    for L in range(1, 9):
        p = mul(power([-1, 1], L), power([-2, 1], L))
        for _ in range(L):
            p = [j * p[j] for j in range(1, len(p))]
        p = scale(p, F(1, factorial(L)))
        ints = [sum(c * F(2 ** (j + r + 1) - 1, j + r + 1) for j, c in enumerate(p)) for r in range(L + 1)]
        check(ints[:L] == [0] * L and ints[-1] == F(factorial(L) ** 2, factorial(2 * L + 1)), "Legendre")
        out.append({"L": L, "first_changed_moment_integral": str(ints[-1])})
    return out


def check_sharpness_family():
    out = []
    for n, beta in ((0, F(1)), (1, F(1)), (0, F(2, 5)), (3, F(1, 2))):
        m = n + 1
        for e in range(1, 31):
            K = 2 ** e
            mom = [F(1) + F(K) ** (2 * m + 1 - 2 * l) + F(factorial(l), K) for l in range(n + 6)]
            p, q = delta_polynomials(mom, n, beta, 3)
            if discriminant_cubic(q) < 0:
                check(gcd_poly(p, q) == [1], "sharpness gcd")
                out.append({"n": n, "beta": str(beta), "first_K": K})
                break
        else:
            raise RuntimeError("sharpness family: no negative discriminant found")
    return out


# ---------------------------------------------------------------------------
# Section C: the all-shift quartic certificate (ported)
# ---------------------------------------------------------------------------

def polynomial_determinant(matrix):
    size = len(matrix)
    dp = {0: [F(1)]}
    for mask in range(1 << size):
        if mask not in dp:
            continue
        row = bin(mask).count("1")
        if row == size:
            continue
        for col in range(size):
            if mask & (1 << col) or trim(matrix[row][col]) == [0]:
                continue
            sign = (-1) ** bin(mask >> (col + 1)).count("1")
            nxt = mask | (1 << col)
            dp[nxt] = add(dp.get(nxt, [F(0)]), scale(mul(dp[mask], matrix[row][col]), sign))
    return trim(dp.get((1 << size) - 1, [F(0)]))


def symbolic_sylvester(p, q):
    dp_, dq = len(p) - 1, len(q) - 1
    width = dp_ + dq
    rows = []
    for coeffs, shifts in ((p, dq), (q, dp_)):
        desc = list(reversed(coeffs))
        for off in range(shifts):
            rows.append([[F(0)]] * off + desc + [[F(0)]] * (width - len(desc) - off))
    return rows


def quartic_discriminant_formula(cf):
    e, d, c, b2, a = cf
    terms = [(256, 3, 0, 0, 0, 3), (-192, 2, 1, 0, 1, 2), (-128, 2, 0, 2, 0, 2), (144, 2, 0, 1, 2, 1),
             (-27, 2, 0, 0, 4, 0), (144, 1, 2, 1, 0, 2), (-6, 1, 2, 0, 2, 1), (-80, 1, 1, 2, 1, 1),
             (18, 1, 1, 1, 3, 0), (16, 1, 0, 4, 0, 1), (-4, 1, 0, 3, 2, 0), (-27, 0, 4, 0, 0, 2),
             (18, 0, 3, 1, 1, 1), (-4, 0, 3, 0, 3, 0), (-4, 0, 2, 3, 0, 1), (1, 0, 2, 2, 2, 0)]
    out = []
    for mult, *pw in terms:
        t = [F(mult)]
        for poly, ex in zip((a, b2, c, d, e), pw):
            t = mul(t, power(poly, ex))
        out.append(t)
    return add(*out)


def quartic_q_in_shift(moments, n):
    """Coefficients of Q_4^(n,beta) as polynomials in s = beta (b = n + s)."""
    q = []
    for j in range(5):
        rf = [F(1)]
        for off in range(3):
            rf = mul(rf, [F(n + j + off), F(1)])
        q.append(scale(rf, F(comb(4, j)) / moments[n + j + 1]))
    return q


def all_shift_certificate():
    moments = quartic_moments(10)
    q = quartic_q_in_shift(moments, 0)
    disc = quartic_discriminant_formula(q)
    res = polynomial_determinant(symbolic_sylvester(q, [scale(q[j], j) for j in range(1, 5)]))
    sdisc, rem = divmod_poly(res, q[-1])
    check(rem == [0] and sdisc == disc, "two symbolic discriminant methods agree")
    check(len(disc) == 19 and all(c < 0 for c in disc), "all 19 coefficients negative")
    neg, mult = primitive(scale(disc, -1))
    fac = [F(1)]
    for cst, ex in ((1, 1), (2, 2), (4, 2), (5, 1)):
        fac = mul(fac, power([F(cst), F(1)], ex))
    cof, r = divmod_poly(neg, fac)
    check(r == [0] and len(cof) == 13 and all(c > 0 and c.denominator == 1 for c in cof), "factorisation")
    hankel = [[moments[1 + i + j] for j in range(3)] for i in range(3)]
    minors = [determinant([row[:s] for row in hankel[:s]]) for s in (1, 2, 3)]
    check(all(v > 0 for v in minors), "Hankel minors")
    # NEW: same measure, starting index n = 1 and n = 2, all shifts beta > 0.
    other = {}
    for n in (1, 2):
        qn = quartic_q_in_shift(moments, n)
        dn = quartic_discriminant_formula(qn)
        signs = sorted(set((c > 0) - (c < 0) for c in dn))
        if n == 1:
            # Proposition prop:allshift for n = 1: degree 18 in beta, all 19 coefficients negative.
            check(len(dn) == 19 and all(c < 0 for c in dn), "n=1: all 19 discriminant coefficients negative")
        other[str(n)] = {"degree_in_beta": len(dn) - 1,
                         "all_coefficients_negative": all(c < 0 for c in dn),
                         "discriminant_at_beta_0_sign": (dn[0] > 0) - (dn[0] < 0),
                         "coefficient_signs_present": signs}
    # Same measure fails at n = 2: Q_4^(2,1) has four (negative) real zeros.
    _, q21 = delta_polynomials(moments, 2, F(1), 4)
    nreal21, sqf21 = real_root_count(q21)
    check(nreal21 == 4 and sqf21 and all(c > 0 for c in q21), "n=2, beta=1: four negative real zeros")
    # Compactness corollary (all shifts, every fixed n): limit quadratic g_b and its discriminant.
    bb = lambda c: [F(c), F(1)]  # b + c
    h = [scale(mul(mul(bb(j), bb(j + 1)), bb(j + 2)), F(comb(4, j))) for j in range(5)]
    check(mul(h[3], bb(6)) == mul(scale(bb(3), F(4)), h[4]), "h_3/h_4 = 4(b+3)/(b+6)")
    check(mul(h[2], mul(bb(5), bb(6))) == mul(scale(mul(bb(2), bb(3)), F(6)), h[4]),
          "h_2/h_4 = 6(b+2)(b+3)/((b+5)(b+6))")
    disc_num = add(scale(mul(power(bb(3), 2), bb(5)), F(16)), scale(mul(mul(bb(2), bb(3)), bb(6)), F(-24)))
    check(disc_num == scale(mul(bb(3), [F(6), F(8), F(1)]), F(-8)),
          "(b+5)(b+6)^2 Disc(g_b) = -8(b+3)(b^2+8b+6)")
    check(add(scale(mul(bb(3), [F(6), F(8), F(1)]), F(10)), scale(mul(bb(5), power(bb(6), 2)), F(-1)))
          == [F(0), F(204), F(93), F(9)], "10(b+3)(b^2+8b+6)-(b+5)(b+6)^2 = 9b^3+93b^2+204b")
    other["2"]["beta_1_real_zero_count"] = nreal21
    return {"measure": QUARTIC_MEASURE, "n": 0, "k": 4,
            "compactness_limit_quadratic": "g_b = z^2 + 4(b+3)/(b+6) z + 6(b+2)(b+3)/((b+5)(b+6))",
            "compactness_disc_identity": "(b+5)(b+6)^2 Disc(g_b) = -8(b+3)(b^2+8b+6) <= -(4/5)(b+5)(b+6)^2",
            "discriminant_degree_in_beta": 18, "all_coefficients_negative": True,
            "negative_discriminant_multiplier": str(mult),
            "factorisation": "-Disc = multiplier*(b+1)(b+2)^2(b+4)^2(b+5)*C(b)",
            "degree_12_cofactor_positive_integer_coefficients_ascending": [int(c) for c in cof],
            "hankel_leading_minors": [str(v) for v in minors],
            "same_measure_other_starting_indices": other}


# ---------------------------------------------------------------------------
# Section D (NEW): persistence for one fixed measure
# ---------------------------------------------------------------------------

def persistence_certificate(kmax):
    mom = cubic_moments(kmax + 4)
    rows = []
    for k in range(1, kmax + 1):
        p, q = delta_polynomials(mom, 0, F(1), k)
        nr = nonreal_count(q)
        g = gcd_degree(p, q)
        check(g == 0, f"gcd(P,Q)=1 at k={k}")
        if k >= 3:
            check(nr >= 2, f"nonreal zeros at k={k}")
        else:
            check(nr == 0, f"real zeros at k={k}")
        rows.append({"k": k, "nonreal_zero_count_Q": nr, "gcd_degree_P_Q": g,
                     "nonreal_zero_count_borghi_P": nonreal_count(borghi_P(mom, 0, k))})
    # Limit function E(w) = sum_j w^j / (j! (b)_j mu_(m+j)), b = m = 1. Jensen polynomials
    # g_d(w) = sum_j C(d,j) w^j / ((1)_j mu_(1+j)); one with nonreal zeros proves E not in LP.
    jensen = []
    first = None
    for d in range(1, 13):
        g = [F(comb(d, j)) / (rising(F(1), j) * mom[1 + j]) for j in range(d + 1)]
        c = nonreal_count(g)
        jensen.append({"d": d, "nonreal_zero_count": c})
        if c > 0 and first is None:
            first = d
    check(first == 6, "first Jensen polynomial of E with nonreal zeros has degree 6")
    g6 = primitive([F(comb(6, j)) / (rising(F(1), j) * mom[1 + j]) for j in range(7)])[0]
    # Order of E is at most 1/2: a_j <= 1/(j!)^2 since mu_l >= 1.
    check(all(m_ >= 1 for m_ in mom), "mu_l >= 1, so E has order <= 1/2")
    return {"measure": CUBIC_MEASURE, "n": 0, "beta": "1", "orders": rows,
            "all_orders_3_to_kmax_have_nonreal_zeros": True, "kmax": kmax,
            "jensen_polynomials_of_limit_E": jensen, "first_nonreal_jensen_degree": first,
            "jensen_g6_primitive_integer_ascending": g6}


# ---------------------------------------------------------------------------
# Section E (NEW): explicit tail-family examples for k = 5..8
# ---------------------------------------------------------------------------

def tail_family_table():
    out = []
    for k in range(5, 9):
        for n in (0, 1):
            for beta in (F(1, 10), F(1), F(5)):
                m = n + 1
                l0 = m + k - 3
                found = None
                for e in range(1, 17):
                    K = 2 ** e
                    mom = [F(1) + F(K) ** (2 * l0 + 1 - 2 * l) + F(factorial(l), K) for l in range(n + k + 3)]
                    p, q = delta_polynomials(mom, n, beta, k)
                    nr = nonreal_count(q)
                    if nr > 0:
                        check(gcd_degree(p, q) == 0, "tail family coprime")
                        found = {"k": k, "n": n, "beta": str(beta), "K": K, "nonreal_zero_count_Q": nr,
                                 "gcd_degree_P_Q": 0,
                                 "borghi_P_nonreal_zero_count": nonreal_count(borghi_P(mom, n, k))}
                        break
                check(found is not None, f"tail family k={k} n={n} beta={beta}")
                out.append(found)
    return out


# ---------------------------------------------------------------------------
# Section F (NEW): natural moment families and the positive theorem
# ---------------------------------------------------------------------------

def natural_families(N):
    fam = {}
    fam["factorial l!"] = [F(factorial(l)) for l in range(N)]
    fam["(2l)!"] = [F(factorial(2 * l)) for l in range(N)]
    fam["Gamma(l+1/2)/sqrt(pi)"] = [rising(F(1, 2), l) for l in range(N)]
    fam["Gamma(l+5)"] = [F(factorial(l + 4)) for l in range(N)]
    fam["1/(l+1) (uniform on [0,1])"] = [F(1, l + 1) for l in range(N)]
    fam["Gamma(3l+1) (superfactorial, nu=3)"] = [F(factorial(3 * l)) for l in range(N)]
    for nu in (F(0), F(1, 4), F(1, 3)):
        fam[f"Bessel K_nu moments, nu={nu}"] = [rising(F(1, 2) + nu, l) * rising(F(1, 2) - nu, l)
                                                / (2 ** l * factorial(l)) for l in range(N)]
    for al, be in ((F(-1, 2), F(1, 2)), (F(1, 2), F(-1, 2)), (F(1, 3), F(2, 3))):
        fam[f"Jacobi weight t^{al}(1-t)^{be} on [0,1]"] = [rising(al + 1, l) / rising(al + be + 2, l)
                                                           for l in range(N)]
    fam["log-normal 2^(l^2)"] = [F(2) ** (l * l) for l in range(N)]
    return fam


def natural_family_tests():
    fams = natural_families(16)
    out = {}
    total = 0
    for name, mom in sorted(fams.items()):
        q_bad, p_bad = [], []
        for k in range(4, 9):
            for n in (0, 1, 2):
                if k <= 8:
                    if distinct_nonreal_count(borghi_P(mom, n, k)) != 0:
                        p_bad.append([k, n])
                for beta in (F(1, 10), F(1, 2), F(1), F(2), F(5)):
                    _, q = delta_polynomials(mom, n, beta, k)
                    total += 1
                    if distinct_nonreal_count(q) != 0:
                        q_bad.append([k, n, str(beta)])
        check(not q_bad, f"natural family {name}: all Q_k real-rooted on the grid")
        out[name] = {"Q_cases_with_nonreal_zeros": q_bad, "borghi_P_cases_with_nonreal_zeros": p_bad}
    # Laguerre identity: for mu_l = Gamma(l+a), P_k is proportional to L_k^(m+a-1)(-z).
    lag = 0
    for a in (F(1), F(1, 2), F(5)):
        for n in (0, 1):
            m = n + 1
            alpha = m + a - 1
            for k in range(1, 9):
                mom_ratio = [rising(F(m) + a, j) for j in range(k + 1)]  # Gamma(m+a+j)/Gamma(m+a)
                p = [F(comb(k, j)) / mom_ratio[j] for j in range(k + 1)]
                lpoly = [rising(alpha + 1 + j, k - j) / factorial(k - j) / factorial(j) for j in range(k + 1)]
                # lpoly_j = C(k+alpha, k-j)/j! * (common factor); compare ratios
                ratios = {p[j] / lpoly[j] for j in range(k + 1)}
                check(len(ratios) == 1, "Laguerre proportionality")
                # exact constant: P_k = k!/Gamma(m+a+k) * L_k^(m+a-1)(-z), i.e. ratio k!/(m+a)_k here
                check(ratios == {F(factorial(k)) / rising(F(m) + a, k)}, "Laguerre constant")
                lag += 1
    # Moment-ratio criterion: delta weight ratios H_j exceed 1 (checked exactly on a grid).
    hmin = None
    for k in range(3, 31):
        for beta in (F(1, 100), F(1), F(10)):
            for n in (0, 3):
                b = n + beta
                h = [comb(k, j) * rising(b + j, k - 1) for j in range(k + 1)]
                for j in range(1, k):
                    r = F(h[j] ** 2, 1) / (h[j - 1] * h[j + 1])
                    check(r > 1, "H_j > 1")
                    hmin = r if hmin is None or r < hmin else hmin
    # Euler moments: the factor-four test fails although Q_k is real-rooted.
    mom = [F(factorial(l)) for l in range(14)]
    _, q = delta_polynomials(mom, 0, F(1), 10)
    ratios = [q[j] ** 2 / (q[j - 1] * q[j + 1]) for j in range(1, 10)]
    check(min(ratios) < 4 and nonreal_count(q) == 0, "Euler: factor-four fails yet real-rooted")
    # Large shift: beta^-(k-1) Q_k -> P_k as beta -> oo. Jacobi weight t^(-1/2)(1-t)^(1/2), n = 0, k = 4,
    # beta = 20 (outside the grid): Q_4 has a genuine nonreal pair (P_4 has nonreal zeros).
    jac = fams["Jacobi weight t^-1/2(1-t)^1/2 on [0,1]"]
    pj, qj = delta_polynomials(jac, 0, F(20), 4)
    qj_int = int_poly(qj)
    check(qj_int == [40425, 371910, 1020096, 1104000, 416000], "Jacobi beta=20 primitive denominator")
    disc_j = quartic_discriminant_formula([[c] for c in qj_int])[0]
    check(disc_j < 0, "Jacobi beta=20: negative discriminant")
    check(nonreal_count(qj) == 2 and gcd_degree(pj, qj) == 0, "Jacobi beta=20: two nonreal zeros, gcd 1")
    check(nonreal_count(borghi_P(jac, 0, 4)) == 2, "Jacobi: P_4 has two nonreal zeros")
    jacobi_large_shift = {"weight": "t^(-1/2)(1-t)^(1/2) on [0,1]", "n": 0, "k": 4, "beta": "20",
                          "primitive_Q_ascending": qj_int, "discriminant": str(disc_j),
                          "nonreal_zero_count_Q": 2, "gcd_degree_P_Q": 0,
                          "nonreal_zero_count_borghi_P": 2}
    return {"families": out, "grid": "k=4..8, n in {0,1,2}, beta in {1/10,1/2,1,2,5}",
            "jacobi_large_shift_example": jacobi_large_shift,
            "number_of_Q_sturm_counts": total, "laguerre_proportionality_cases": lag,
            "min_delta_weight_ratio_H_j_on_grid": r8(hmin),
            "euler_k10_min_factor_four_ratio": r8(min(ratios)),
            "euler_k10_nonreal_zero_count": 0}


# ---------------------------------------------------------------------------
# Section G: numerical illustrations (mpmath; not part of any proof)
# ---------------------------------------------------------------------------

def numerical_sections():
    try:
        import mpmath as mp
    except ImportError:
        return {"status": "skipped: mpmath not installed"}
    mp.mp.dps = 60

    def cr(z):
        return [r8(mp.re(z)), r8(mp.im(z))]

    out = {}
    # G1: cubic example poles and residues.
    mom = cubic_moments(9)
    p, q = delta_polynomials(mom, 0, F(1), 3)
    pi, _ = primitive(p)
    qi, _ = primitive(q)
    roots = mp.polyroots([mp.mpf(c) for c in reversed(qi)], maxsteps=200, extraprec=200)
    dq = [j * qi[j] for j in range(1, len(qi))]
    res = []
    for r in sorted(roots, key=lambda z: (float(mp.re(z)), float(mp.im(z)))):
        val = mp.polyval([mp.mpf(c) for c in reversed(pi)], r) / (
            1000 * r * mp.polyval([mp.mpf(c) for c in reversed(dq)], r))
        res.append({"pole": cr(r), "residue": cr(val)})
    res0 = mp.mpf(pi[0]) / (1000 * qi[0])

    def f_true(z):
        return 1 / (z + 1) + 50 / (z + mp.mpf(1) / 50) + mp.quad(lambda t: mp.e ** (-t) / (z + t), [0, mp.inf]) / 1000

    def d3(z):
        return mp.polyval([mp.mpf(c) for c in reversed(pi)], z) / (1000 * z * mp.polyval([mp.mpf(c) for c in reversed(qi)], z))

    errs = {}
    for label, z in (("z=10", mp.mpf(10)), ("z=1", mp.mpf(1)), ("z=0.1", mp.mpf("0.1")), ("z=i", mp.mpc(0, 1)),
                     ("z=-0.5+0.3i", mp.mpc(-0.5, 0.3))):
        errs[label] = r8(abs(d3(z) - f_true(z)) / abs(f_true(z)))
    out["numerical_cubic_example_poles_residues"] = {"poles_and_residues": res, "residue_at_z0": r8(res0),
                                                     "relative_errors": errs}
    # G2: persistence: scaled zeros k^2 z near the limit zero w0 of E.
    a = [1 / (mp.factorial(j) ** 2 * (1 + mp.mpf(50) ** (-j) + mp.factorial(j + 1) / 1000)) for j in range(90)]
    w0 = mp.findroot(lambda w: mp.fsum(a[j] * w ** j for j in range(90)), mp.mpc(-9.1, 9.13))
    scaled = {}
    momk = cubic_moments(40)
    for k in (10, 20, 30):
        _, qk = delta_polynomials(momk, 0, F(1), k)
        qk_int, _ = primitive(qk)
        rts = mp.polyroots([mp.mpf(c) for c in reversed(qk_int)], maxsteps=800, extraprec=3000)
        best = min((r * k * k for r in rts), key=lambda w: abs(w - w0))
        scaled[str(k)] = cr(best)
    out["numerical_persistence_scaled_zeros"] = {"limit_zero_w0_of_E": cr(w0), "nearest_k2_z_by_order": scaled}
    # G3: co-occurrence of nonreal poles and accuracy loss for f = 100/(z+2) + int e^-t/(z+t) dt at z = 1.
    momx = [F(100 * 2 ** l + factorial(l)) for l in range(60)]
    ftrue = mp.mpf(100) / 3 + mp.e * mp.e1(1)
    rows = []
    for k in (2, 4, 8, 12, 16, 20, 30, 40, 50):
        pk, qk = delta_polynomials(momx, 0, F(1), k)
        val = F(evaluate(pk, F(1))) / evaluate(qk, F(1))
        err = abs(mp.mpf(val.numerator) / val.denominator - ftrue) / ftrue
        row = {"k": k, "relative_error_at_z1": r8(err), "nonreal_zero_count_Q_exact": nonreal_count(qk),
               "gcd_degree_P_Q_exact": gcd_degree(pk, qk)}
        rows.append(row)
    momE = [F(factorial(l)) for l in range(60)]
    ftrueE = mp.e * mp.e1(1)
    for k in (8, 30):
        pk, qk = delta_polynomials(momE, 0, F(1), k)
        val = F(evaluate(pk, F(1))) / evaluate(qk, F(1))
        rows.append({"k": k, "series": "Euler only", "relative_error_at_z1":
                     r8(abs(mp.mpf(val.numerator) / val.denominator - ftrueE) / ftrueE)})
    out["numerical_accuracy_cooccurrence"] = {"function": "100/(z+2) + int_0^inf exp(-t)/(z+t) dt",
                                              "moments": "100*2^l + l!", "n": 0, "beta": "1", "rows": rows,
                                              "note": "co-occurrence only; no causal claim"}
    return out


# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--kmax", type=int, default=30, help="highest order in the persistence certificate")
    args = parser.parse_args()
    report = {
        "A_symbolic_identities": check_cubic_identities(),
        "A_deltaQ_tailratio_rational_form": check_deltaQ_and_tailratio(),
        "B_cubic_example": exact_cubic_example(),
        "B_quartic_example": exact_quartic_example(),
        "B_legendre_perturbation": check_legendre(),
        "B_cubic_sharpness_family": check_sharpness_family(),
        "C_all_shift_quartic": all_shift_certificate(),
        "D_persistence_fixed_measure": persistence_certificate(args.kmax),
        "E_tail_family_explicit_examples": tail_family_table(),
        "F_natural_families_and_positive_theorem": natural_family_tests(),
    }
    report.update(numerical_sections())
    report["status"] = "all exact checks passed"
    report["scope"] = ("Exact sections certify the finite statements in the paper; all-parameter statements are "
                       "proved in the paper. numerical_* sections are floating-point illustrations only.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=1, sort_keys=True) + "\n"
    (args.output_dir / "verification.json").write_text(text, encoding="utf-8")
    print("all exact checks passed")
    print("persistence orders certified: 3..%d" % args.kmax)
    return 0


if __name__ == "__main__":
    sys.exit(main())
