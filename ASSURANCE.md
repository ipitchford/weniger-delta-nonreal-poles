# Assurance

This file says what has been checked and what has not. Each row is a separate kind of assurance, and no row implies another.

| Dimension | State | Evidence |
|---|---|---|
| Availability | Public GitHub release and Zenodo DOI 10.5281/zenodo.23586255 | Release assets and [SHA256SUMS](SHA256SUMS) |
| Written proofs | Complete proofs in the manuscript for every result marked "proof" in Table 2 | [paper/paper.pdf](paper/paper.pdf). The limit and Hurwitz arguments (Section 6), the Rouché arguments (Theorem 3.2, Corollary 5.2), the perturbation and smoothing arguments (Lemma 4.5, Corollary 4.6) and the positive theorem (Section 7) are proved in prose only. |
| Exact certificates | Exact rational arithmetic throughout sections A–F | [code/verify_all.py](code/verify_all.py) uses Python integers and `fractions.Fraction` only: primitive pseudo-remainder Sturm sequences, exact gcds, Sylvester determinants and symbolic identities over Q[b], Q[β] and Q[b, x, y]. Every check is an explicit `RuntimeError`, so `python -O` disables nothing. |
| What the code checks | Stated in Appendix B of the manuscript | The discriminant identity (12) over Q[b, x, y]; identity (3) for n ≤ 3, k ≤ 8; identity (8) for 3 ≤ k ≤ 12 and the factorisation (k − 4)(4k − 5); the rational form (4) against definition (2) at 81 rational points; the explicit examples; the Legendre integrals; the cubic sharpness family; Proposition 5.1 for n = 0 and n = 1 (both asserted), the four real zeros at n = 2, and the identities for the limit quadratic of Corollary 5.2; Theorem 6.3(1) and the Jensen polynomials g_1, …, g_12; Table 3; the natural-family grid, the Laguerre identity with its constant, and the β = 20 Jacobi example. Finite checks supplement the proofs; they do not replace the all-parameter arguments. |
| Numerical results | Floating point (mpmath, 60 digits), labelled as such | `numerical_*` sections of [results/verification.json](results/verification.json); Section 8 and Figure 1. |
| Internal replay | Passed | Fresh replays of the final package under `python -B -O` reproduced `results/verification.json` byte for byte (SHA-256 `c0cb9b0c…3cb3bab84`). Details in [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). |
| Semantic negative controls | Passed | [negative_controls.py](negative_controls.py) runs an unmodified baseline, then checks that each of the following is rejected: a wrong cubic threshold constant; a shift above the threshold claimed safe; a perturbed counterexample measure; flipped Sturm counts in Table 4 and Table 3; three tampered result values (including the β = 20 Jacobi certificate); a gcd routine that reports every pair coprime (caught by an injected common factor); three code mutations (rising-factorial length in the denominator, a cubic discriminant coefficient, un-negated Sturm remainders); and a tampered file. |
| Implementation diversity | Partial | The first-round reviewers recomputed the printed integers, the corner factorisations, b_c, the tail-ratio bound, the degree-18 discriminant and the Hankel minors in sympy and mpmath. Both adversarial audits wrote their own exact SymPy code without importing package code: on the first draft, 48 general-order witnesses, 48 safe cubic quadruples, 42 single-inequality isolation pairs and a symbolic resultant for the all-shift quartic; on the revised version, the genuine-pole counts for k = 1–30, the Jensen counts to degree 12, the all-shift quartic at n = 0, 1, 2, all 24 tail-family cases, all 975 grid cases plus 180 further cases, and the algebra of Corollary 5.2. The re-review ran its own SymPy spot checks. This is diversity inside one producer workflow, not independent reimplementation. |
| Internal AI-assisted review | Done | First-round multi-perspective review and an AI-assisted adversarial audit by a different model (configured as OpenAI gpt-6-astra), both of the 19 October draft; a verification re-review of the revised version; and a second adversarial audit of the revised version by the same external model, followed by a final editor pass. Summary in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). This is not peer review. |
| Independent rerun or reimplementation by an unaffiliated party | Not established | |
| Specialist (human) review | Not established | |
| Editorial peer review | Not established | |
| Formal verification | Not attempted | No proof assistant was used. |
| Source fidelity | Published article and author preprint | Borghi's question and the §3.3 sentences are quoted from the published article; the re-review and the second audit rechecked every quoted fragment against freshly fetched primary text. The formulas, absent from the HTML, were checked against the author preprint (equations 26, 29, 30, 35, 36). |
| Environment | Standard library plus mpmath 1.3.0; deterministic JSON | [ENVIRONMENT.txt](ENVIRONMENT.txt), [requirements.txt](requirements.txt) |
| Novelty | Bounded search only | [RESEARCH_GATES.json](RESEARCH_GATES.json), [SOURCES.md](SOURCES.md) |
| Citations | Metadata checked | [CITATION_AUDIT.md](CITATION_AUDIT.md) |

## Scope and limits

* All results are at fixed finite order, and Theorem 4.1 is a worst case over all admissible measures: the witness measure depends on (k, n, β). Section 6 is the only fixed-measure result, and it gives genuine poles only for 3 ≤ k ≤ 30. It rules out an argument that requires the unreduced denominator to be real-rooted at all large orders; it does not show that genuine nonreal poles persist.
* The threshold k_0 in Theorem 6.3(2) is not effective, and Proposition 4.3 does not exclude cancellation at large k. Conjecture 6.4 is open.
* The explicit all-shift measure of Proposition 5.1 is certified for n = 0 and n = 1 only; Corollary 5.2 gives existence for each fixed n without explicit constants.
* The natural-family grid (Section 7) certifies real-rootedness at the grid points only (β ≤ 5). For the K_ν weights it is evidence, not proof; the non-integer Jacobi weights fail at larger shifts (Remark 7.4).
* Hutchinson (1923) and Kurtz (1992) are cited without theorem locators; the strict criterion is proved in the manuscript (Lemma 3.1).
* Byte identity of the `numerical_*` sections is expected only with mpmath 1.3.0; they are floating point.

Byte-identical replay shows that the shipped output is exactly what the shipped program produces. It does not by itself show that the program implements the manuscript's definitions. That rests on the check of the rational form against the finite-difference definition (81 points), the negative controls, and the independent recomputations in the internal reviews.
