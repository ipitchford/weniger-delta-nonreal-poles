#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Author: Anonymous
"""Semantic negative controls for this package.

A replay that cannot fail is not evidence.  Each control below copies code/ and results/
into a fresh temporary directory, corrupts one claim, one line of code or one output, and
requires an exact checker to REJECT it.  An unmodified baseline must pass first.  Nothing
is written beside this script or anywhere in the package directory.

The checker (run as a child process `python -B -O negative_controls.py --child PKG CLAIMS`)
imports the package's own exact code (code/verify_all.py) from the copy and requires:
  * the exact sections A, B, C, E and F recomputed from the code to equal the corresponding
    sections of results/verification.json; section D (persistence) recomputed to order 14 to
    equal the first 14 shipped rows, the Jensen data and the g_6 coefficients;
  * the cubic threshold claim: at b just below the claimed constant both internal factor-four
    inequalities of Q_3 hold (Theorem 3.2 then gives universal safety), and just above it one
    fails;
  * every shift claimed safe at order three: both inequalities hold, and the cubic sharpness
    family (Appendix A) has a positive discriminant for K = 2, 4, ..., 2^30;
  * every claimed counterexample measure (atoms plus c e^{-t} dt): the claimed exact number of
    nonreal zeros of Q_k, at least two, and gcd(P, Q) = 1;
  * the claimed Table 4 counts (k = 3..14) and Table 3 entries recomputed by exact Sturm counts;
  * a gcd regression: multiplying the cubic example's P and Q by z^2 + z + 1 must give a gcd of
    degree 2 (so a gcd routine that reports every pair coprime is caught).

Controls:
  1. baseline: the unmodified package passes the checker;
  2. wrong threshold constant: b_c claimed to be the positive root of b^2 + 3b - 6
     (= b_c + 1, the other corner) instead of b^2 + 5b - 2 -> rejected;
  3. wrong safe-side shift: n = 0, beta = 2/5 (> b_c) claimed universally safe at k = 3 -> rejected;
  4. perturbed counterexample: the cubic measure with the atom 50 delta_(1/50) replaced by
     2 delta_(1/2), still claimed to give nonreal zeros at k = 3 -> rejected (Q_3 is real-rooted);
  5. flipped Sturm counts: Table 4 at k = 8 claimed to have 2 nonreal zeros (true: 4); Table 3 at
     k = 5, n = 0, beta = 1 claimed to need only K = 2^9 (true: 2^10) -> each rejected;
  6. tampered outputs in results/verification.json: one digit of the cubic discriminant; the
     persistence count at k = 8; the nonreal count of the Jacobi beta = 20 example of Remark 7.4
     -> each rejected;
  7. gcd != 1 injection: gcd_degree mutated to report every pair coprime -> rejected by the
     injected common factor;
  8. mutated identities in code/verify_all.py: the rising factorial (b+j)_(k-1) in the
     denominator replaced by (b+j)_k; the coefficient 18 of the cubic discriminant replaced by 16;
     the Sturm sequence built from un-negated remainders -> each rejected;
  9. a tampered file is detected by SHA256SUMS (shipped sums if present, else sums of the copy).

Usage:  python -B -O negative_controls.py      Exit status 0 iff every control behaves.
Standard library only.  Runtime a few seconds to about half a minute.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.abspath(__file__))

CUBIC_ATOMS = [["1", "1"], ["50", "1/50"]]          # [weight, location]: delta_1 + 50 delta_(1/50)
TABLE4 = {3: (2, 2), 4: (2, 2), 5: (2, 2), 6: (2, 2), 7: (2, 2), 8: (4, 4), 9: (4, 6), 10: (4, 8),
          11: (6, 8), 12: (6, 8), 13: (6, 10), 14: (6, 10)}
TABLE3 = [  # (k, n, beta, K, nonreal count of Q_k)
    (5, 0, "1/10", 2 ** 11, 2), (5, 0, "1", 2 ** 10, 2), (5, 0, "5", 2, 2),
    (6, 1, "1", 2, 2), (7, 0, "5", 2, 4), (8, 1, "5", 2, 6),
]
BASE_CLAIMS = {
    "threshold_quadratic_ascending": [-2, 5, 1],       # b_c = positive root of b^2 + 5b - 2
    "safe_cubic_shifts": [[0, "3/10"], [0, "1/3"]],
    "counterexamples": [
        {"label": "cubic example, Section 5", "atoms": CUBIC_ATOMS, "density": "1/1000",
         "n": 0, "beta": "1", "k": 3, "nonreal": 2},
        {"label": "quartic example, Section 5", "atoms": [["1", "1"], ["100000", "1/100"]], "density": "1/1000",
         "n": 0, "beta": "1/3", "k": 4, "nonreal": 2},
    ],
    "table4": {str(k): list(v) for k, v in TABLE4.items()},
    "table3": [list(r) for r in TABLE3],
}
COMPARED_SECTIONS = ["A_symbolic_identities", "A_deltaQ_tailratio_rational_form", "B_cubic_example",
                     "B_quartic_example", "B_legendre_perturbation", "B_cubic_sharpness_family",
                     "C_all_shift_quartic", "E_tail_family_explicit_examples",
                     "F_natural_families_and_positive_theorem"]
D_ORDERS = 14


# ------------------------------------------------------------------ child: the exact checker
def child(pkg, claims_path):
    sys.path.insert(0, os.path.join(pkg, "code"))
    from fractions import Fraction as F
    from math import comb, factorial
    import verify_all as V

    errs = []
    with open(claims_path) as fh:
        claims = json.load(fh)
    with open(os.path.join(pkg, "results", "verification.json")) as fh:
        shipped = json.load(fh)

    def norm(x):
        return json.loads(json.dumps(x, sort_keys=True))

    def guarded(label, fn):
        try:
            return fn()
        except RuntimeError as exc:
            errs.append("REJECT %s: %s" % (label, exc))
            return None

    if shipped.get("status") != "all exact checks passed":
        errs.append("MISMATCH: status")
    sections = {
        "A_symbolic_identities": V.check_cubic_identities,
        "A_deltaQ_tailratio_rational_form": V.check_deltaQ_and_tailratio,
        "B_cubic_example": V.exact_cubic_example,
        "B_quartic_example": V.exact_quartic_example,
        "B_legendre_perturbation": V.check_legendre,
        "B_cubic_sharpness_family": V.check_sharpness_family,
        "C_all_shift_quartic": V.all_shift_certificate,
        "E_tail_family_explicit_examples": V.tail_family_table,
        "F_natural_families_and_positive_theorem": V.natural_family_tests,
    }
    for name in COMPARED_SECTIONS:
        got = guarded(name, sections[name])
        if got is not None and norm(got) != shipped.get(name):
            errs.append("MISMATCH: " + name)
    d = guarded("D_persistence", lambda: V.persistence_certificate(D_ORDERS))
    sd = shipped.get("D_persistence_fixed_measure", {})
    if d is not None:
        d = norm(d)
        for key in ("jensen_polynomials_of_limit_E", "first_nonreal_jensen_degree",
                    "jensen_g6_primitive_integer_ascending", "measure", "n", "beta"):
            if d.get(key) != sd.get(key):
                errs.append("MISMATCH: D_persistence_fixed_measure." + key)
        if d["orders"] != sd.get("orders", [])[:D_ORDERS]:
            errs.append("MISMATCH: D_persistence_fixed_measure.orders[:%d]" % D_ORDERS)

    # Cubic weights h_j = C(3,j)(b+j)_2; the two internal factor-four discriminants.
    def cubic_conditions(b):
        h = [comb(3, j) * V.rising(b + j, 2) for j in range(4)]
        return h[1] ** 2 - 4 * h[0] * h[2], h[2] ** 2 - 4 * h[1] * h[3]

    # Threshold: bisect the claimed quadratic's positive root in exact rationals.
    qa = [F(c) for c in claims["threshold_quadratic_ascending"]]
    lo, hi = F(0), F(10)
    if not (V.evaluate(qa, lo) < 0 < V.evaluate(qa, hi)):
        errs.append("REJECT threshold: claimed quadratic has no sign change on (0, 10)")
    else:
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if V.evaluate(qa, mid) < 0 else (lo, mid)
        eps = F(1, 10 ** 9)
        below, above = cubic_conditions(lo - eps), cubic_conditions(hi + eps)
        if not all(c >= 0 for c in below):
            errs.append("REJECT threshold: a factor-four inequality fails just below the claimed constant")
        if all(c >= 0 for c in above):
            errs.append("REJECT threshold: both factor-four inequalities hold just above the claimed constant")

    for n, beta in claims["safe_cubic_shifts"]:
        beta = F(beta)
        if not all(c >= 0 for c in cubic_conditions(n + beta)):
            errs.append("REJECT safe shift n=%d beta=%s: factor-four inequality fails" % (n, beta))
        m = n + 1
        for e in range(1, 31):
            K = 2 ** e
            mom = [F(1) + F(K) ** (2 * m + 1 - 2 * l) + F(factorial(l), K) for l in range(n + 6)]
            _, q = V.delta_polynomials(mom, n, beta, 3)
            if V.discriminant_cubic(q) <= 0:
                errs.append("REJECT safe shift n=%d beta=%s: sharpness family K=%d is not real-rooted" % (n, beta, K))
                break

    for ce in claims["counterexamples"]:
        n, beta, k = ce["n"], F(ce["beta"]), ce["k"]
        dens = F(ce["density"])
        mom = [sum(F(w) * F(a) ** l for w, a in ce["atoms"]) + dens * factorial(l) for l in range(n + k + 3)]
        p, q = V.delta_polynomials(mom, n, beta, k)
        nr = guarded(ce["label"], lambda: V.nonreal_count(q))
        if nr is not None and (nr != ce["nonreal"] or nr < 2):
            errs.append("REJECT %s: Q_%d has %d nonreal zeros, claimed %d" % (ce["label"], k, nr, ce["nonreal"]))
        if V.gcd_degree(p, q) != 0:
            errs.append("REJECT %s: P and Q are not coprime" % ce["label"])

    mom = V.cubic_moments(D_ORDERS + 4)
    for k, (cq, cp) in claims["table4"].items():
        k = int(k)
        _, q = V.delta_polynomials(mom, 0, F(1), k)
        got = guarded("table4", lambda: (V.nonreal_count(q), V.nonreal_count(V.borghi_P(mom, 0, k))))
        if got is not None and got != (cq, cp):
            errs.append("REJECT Table 4 at k=%d: exact counts %s, claimed %s" % (k, got, (cq, cp)))

    for k, n, beta, K, count in claims["table3"]:
        beta = F(beta)
        l0 = n + 1 + k - 3

        def tail_count(KK):
            mm = [F(1) + F(KK) ** (2 * l0 + 1 - 2 * l) + F(factorial(l), KK) for l in range(n + k + 3)]
            return V.nonreal_count(V.delta_polynomials(mm, n, beta, k)[1])
        got = guarded("table3", lambda: (tail_count(K), K == 2 or tail_count(K // 2) == 0))
        if got is not None and got != (count, True):
            errs.append("REJECT Table 3 at k=%d n=%d beta=%s: claimed smallest K=%d with %d nonreal zeros"
                        % (k, n, beta, K, count))

    # gcd regression: inject a known common factor.
    ver = shipped.get("B_cubic_example", {})
    pi, qi = ver.get("P_integer_ascending"), ver.get("Q_integer_ascending")
    if pi and qi:
        common = [1, 1, 1]
        if V.gcd_degree(V.mul(pi, common), V.mul(qi, common)) != 2:
            errs.append("REJECT gcd regression: injected common factor z^2+z+1 not detected")
        if V.gcd_degree(pi, qi) != 0:
            errs.append("REJECT gcd regression: shipped cubic pair not coprime")
    else:
        errs.append("MISMATCH: cubic example coefficients missing")

    for e in errs:
        print(e, file=sys.stderr)
    print("CHECKER " + ("PASS" if not errs else "FAIL"))
    return 0 if not errs else 1


# ------------------------------------------------------------------ parent
def fresh_copy():
    tmp = tempfile.TemporaryDirectory(prefix="bw-negctl-")
    dst = os.path.join(tmp.name, "pkg")
    os.makedirs(dst)
    for sub in ("code", "results"):
        shutil.copytree(os.path.join(ROOT, sub), os.path.join(dst, sub),
                        ignore=shutil.ignore_patterns("__pycache__"))
    return tmp, dst


def run_checker(pkg, claims):
    path = os.path.join(os.path.dirname(pkg), "claims.json")
    with open(path, "w") as fh:
        json.dump(claims, fh)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, "-B", "-O", os.path.abspath(__file__), "--child", pkg, path],
                       capture_output=True, text=True, env=env, cwd=os.path.dirname(pkg))
    return r.returncode, (r.stdout + r.stderr).strip()


def mutate(pkg, old, new):
    path = os.path.join(pkg, "code", "verify_all.py")
    with open(path) as fh:
        src = fh.read()
    if src.count(old) != 1:
        raise SystemExit("control setup error: mutation anchor not unique: %r" % old)
    with open(path, "w") as fh:
        fh.write(src.replace(old, new))


def edit_results(pkg, fn):
    path = os.path.join(pkg, "results", "verification.json")
    with open(path) as fh:
        data = json.load(fh)
    fn(data)
    with open(path, "w") as fh:
        fh.write(json.dumps(data, indent=1, sort_keys=True) + "\n")


def claims_with(**changes):
    c = json.loads(json.dumps(BASE_CLAIMS))
    for key, fn in changes.items():
        fn(c)
    return c


def sums_of(base, names):
    return {n: hashlib.sha256(open(os.path.join(base, n), "rb").read()).hexdigest() for n in names}


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "--child":
        return child(sys.argv[2], sys.argv[3])

    results = []

    def control(label, expect_pass, setup=None, claims=None):
        tmp, pkg = fresh_copy()
        try:
            if setup:
                setup(pkg)
            rc, out = run_checker(pkg, claims if claims is not None else BASE_CLAIMS)
        finally:
            tmp.cleanup()
        ok = (rc == 0) if expect_pass else (rc != 0)
        lines = [l for l in out.splitlines() if not l.startswith("CHECKER")]
        reason = lines[0] if lines else (out.splitlines()[0] if out else "")
        results.append((label, ok))
        print("%-4s %s%s" % ("OK" if ok else "BAD", label, "" if expect_pass else "  [" + reason[:110] + "]"))
        if expect_pass and not ok:
            print(out)

    control("baseline: unmodified package passes", True)
    control("wrong threshold constant (root of b^2+3b-6) rejected", False,
            claims=claims_with(t=lambda c: c.__setitem__("threshold_quadratic_ascending", [-6, 3, 1])))
    control("wrong safe-side shift (n=0, beta=2/5 claimed safe at k=3) rejected", False,
            claims=claims_with(s=lambda c: c["safe_cubic_shifts"].append([0, "2/5"])))
    control("perturbed counterexample measure (2 delta_(1/2)) rejected", False,
            claims=claims_with(p=lambda c: c["counterexamples"][0].__setitem__("atoms", [["1", "1"], ["2", "1/2"]])))
    control("flipped Sturm count, Table 4 k=8 (claimed 2) rejected", False,
            claims=claims_with(f=lambda c: c["table4"].__setitem__("8", [2, 4])))
    control("flipped Table 3 entry, k=5 n=0 beta=1 (claimed K=2^9) rejected", False,
            claims=claims_with(f=lambda c: c["table3"][1].__setitem__(3, 2 ** 9)))

    def tamper_disc(d):
        s = d["B_cubic_example"]["discriminant_Q"]
        d["B_cubic_example"]["discriminant_Q"] = s[:-2] + ("1" if s[-2] != "1" else "2") + s[-1]
    control("tampered output: cubic discriminant digit rejected", False,
            setup=lambda pkg: edit_results(pkg, tamper_disc))

    def tamper_count(d):
        for r in d["D_persistence_fixed_measure"]["orders"]:
            if r["k"] == 8:
                r["nonreal_zero_count_Q"] = 2
    control("tampered output: persistence count at k=8 rejected", False,
            setup=lambda pkg: edit_results(pkg, tamper_count))
    def tamper_jacobi(d):
        d["F_natural_families_and_positive_theorem"]["jacobi_large_shift_example"]["nonreal_zero_count_Q"] = 0
    control("tampered output: Jacobi beta=20 genuine-pole count (Remark 7.4) rejected", False,
            setup=lambda pkg: edit_results(pkg, tamper_jacobi))
    control("gcd != 1 injection: gcd routine reporting every pair coprime rejected", False,
            setup=lambda pkg: mutate(pkg, "        a, b = b, (iprim(r) if any(r) else [0])\n    return len(a) - 1",
                                     "        a, b = b, (iprim(r) if any(r) else [0])\n    return 0"))
    control("mutated identity: (b+j)_(k-1) -> (b+j)_k in the denominator rejected", False,
            setup=lambda pkg: mutate(pkg, "q = [F(comb(k, j)) * rising(F(n) + beta + j, k - 1) / moments[n + j + 1]",
                                     "q = [F(comb(k, j)) * rising(F(n) + beta + j, k) / moments[n + j + 1]"))
    control("mutated identity: cubic discriminant coefficient 18 -> 16 rejected", False,
            setup=lambda pkg: mutate(pkg, "27 * a * a * d * d + 18 * a * b2 * c * d",
                                     "27 * a * a * d * d + 16 * a * b2 * c * d"))
    control("mutated Sturm sequence: remainders not negated rejected", False,
            setup=lambda pkg: mutate(pkg, "r = [-c for c in iprem(seq[-2], seq[-1])]",
                                     "r = [c for c in iprem(seq[-2], seq[-1])]"))

    # Tampered file versus SHA256SUMS (no child needed).
    tmp, pkg = fresh_copy()
    try:
        names = sorted(os.path.relpath(os.path.join(dp, f), pkg).replace(os.sep, "/")
                       for dp, _, fs in os.walk(pkg) for f in fs)
        shipped = os.path.join(ROOT, "SHA256SUMS")
        if os.path.exists(shipped):
            ref = {}
            for line in open(shipped).read().splitlines():
                h, n = line.split(None, 1)
                if n in names:
                    ref[n] = h
            baseline_ok = ref == sums_of(pkg, names)
        else:
            ref = sums_of(pkg, names)
            baseline_ok = True
        with open(os.path.join(pkg, "code", "verify_all.py"), "a") as fh:
            fh.write("# tampered\n")
        detected = ref != sums_of(pkg, names)
    finally:
        tmp.cleanup()
    ok = baseline_ok and detected
    results.append(("tampered file detected by SHA256SUMS", ok))
    print("%-4s %s" % ("OK" if ok else "BAD", "tampered file detected by SHA256SUMS"
                       + ("" if baseline_ok else " (baseline hashes do not match SHA256SUMS)")))

    bad = [label for label, ok in results if not ok]
    print("ALL NEGATIVE CONTROLS BEHAVED" if not bad else "FAILED CONTROLS: %s" % bad)
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
