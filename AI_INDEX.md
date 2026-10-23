# Agent navigation

## Identity

- Release slug: `weniger-delta-nonreal-poles`. Version `1.0.0-candidate`, tag `v1.0.0-candidate`. DOI 10.5281/zenodo.23586255; concept DOI 10.5281/zenodo.23586254.
- Repository: https://github.com/ipitchford/weniger-delta-nonreal-poles. Release page: https://evidencepress.org/releases/weniger-delta-nonreal-poles/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](paper/paper.pdf) (15 pages) and [LaTeX source](paper/paper.tex). The figure is [paper/figures/cubic_nonreal_poles.pdf](paper/figures/cubic_nonreal_poles.pdf) (floating-point coordinates; illustration only).
- Machine-readable claims: [CLAIMS.json](CLAIMS.json). Overview for humans: [README.md](README.md).
- Assurance limits: [ASSURANCE.md](ASSURANCE.md). Internal review summary: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
- Provenance: [PROVENANCE.md](PROVENANCE.md). Sources: [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [RESEARCH_GATES.json](RESEARCH_GATES.json).
- Citation: [CITATION.cff](CITATION.cff) and [Zenodo metadata](.zenodo.json). Licences: [LICENSES.md](LICENSES.md), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md).

## Scientific map

Object: a positive Borel measure μ on [0, ∞) with infinite support and finite moments μ_ℓ ("admissible"); the Stieltjes function f(z) = ∫(z + t)⁻¹ dμ(t); Weniger's delta transformation δ_k^{(n)}(β; z) of its asymptotic series, with starting index n ≥ 0, order k ≥ 1 and shift β > 0. With b = n + β and m = n + 1 the denominator is Q_k(z) = Σ_j C(k, j)(b + j)_{k−1} z^j/μ_{m+j}. Borghi's auxiliary polynomial is P_k(z) = Σ_j C(k, j) z^j/μ_{m+j}. Source: Borghi, Mathematics 14 (2026) 376, doi:10.3390/math14020376, §3.1 and §3.3. The manuscript's β is the published article's γ.

- Borghi's question (all denominator zeros on the negative real axis) is **false** in general, already at n = 0, β = 1, k = 3. Each claim below gives its status and its location in [paper/paper.pdf](paper/paper.pdf).
  - Classification (Table 1): orders 1 and 2 always safe; order 3 safe iff b ≤ b_c = (√33 − 5)/2 ≈ 0.3723; order ≥ 4 never safe. Unsafe means some admissible measure (moments ~ c·ℓ!, or a smooth positive density) gives a genuine nonreal pair of poles. Theorem 4.1, Corollary 4.6, **proved**. The measure depends on (k, n, β).
  - Fixed weights h_j: real-rooted for every admissible measure iff h_j² ≥ 4h_{j−1}h_{j+1}. Theorem 3.2, **proved** (sufficiency is Hutchinson–Kurtz).
  - Genuine poles: accuracy through order, deg gcd(P, Q) < ⌈k/2⌉, and a one-moment perturbation. Lemma 4.2, Proposition 4.3, Lemmas 4.4–4.5, **proved**.
- Explicit examples (Section 5), **exact certificates**: the cubic example (n = 0, β = 1, μ = δ_1 + 50δ_{1/50} + 10⁻³e^{−t}dt); the quartic example (n = 0, β = 1/3); one measure λ defeating every shift at k = 4 for n = 0 and n = 1 (Proposition 5.1; it fails at n = 2); the tail family at k = 5–8 (Table 3). Every counterexample has non-real-rooted P_k (Remark 5.3).
- For every fixed n, one measure defeats every shift at k = 4 (Corollary 5.2, **proved**, no explicit constants).
- One fixed measure (Section 6): Q_k(w/k²)/(b)_{k−1} → E_μ(w) (Theorem 6.1, **proved**); nonreal denominator zeros at every k ≥ k_0 for the cubic measure (Theorem 6.3(2), **proved**, k_0 not effective); genuine nonreal poles at every 3 ≤ k ≤ 30 (Theorem 6.3(1), Table 4, **exact certificate**); every k ≥ 3 (Conjecture 6.4, **open**).
- Positive side (Section 7): real-rootedness at every order if P_k is real-rooted, if (1/μ_{m+j}) is a multiplier sequence, or if μ_{r−1}μ_{r+1} ≥ 4μ_r² (Theorem 7.2, **proved**). Covers Γ(ℓ + a), Γ(νℓ + c) with integer ν, Beta weights with integer exponent, log-normal with σ² ≥ 2 log 2 (Corollary 7.3, **proved**), extending to every n, β and k the cases Borghi proved at n = 0, γ = 1 (§§4.2, 4.3, 4.6). Q_k real-rooted for every shift iff P_k is (Remark 7.4, **proved**); a non-integer Jacobi weight at β = 20 has genuine nonreal poles (**exact certificate**). A finite grid with β ≤ 5: 975 Sturm counts over 15 families (**exact certificate on the grid**; evidence only for K_ν).
- Practice (Section 8): pole locations, complex residues (so δ_3 is not a Stieltjes function), relative errors and an accuracy co-occurrence table, **numerical**; an exact Sturm-count diagnostic.

Do not read:
- Theorem 4.1 as one measure failing at every order (the measure depends on (k, n, β); Section 6 is the fixed-measure result);
- Theorem 6.3 as proving genuine poles beyond k = 30 (beyond that only nonreal denominator zeros are proved);
- Proposition 5.1 as covering n ≥ 2 (Corollary 5.2 gives existence, not an explicit measure);
- any result as a statement about convergence or divergence of the delta approximants;
- the natural-family grid as a proof for the K_ν family, or as holding for shifts above 5;
- the factor-four test as a practical diagnostic (it fails on the Euler moments at k = 10 although Q_10 is real-rooted);
- internal AI-assisted review as peer review.

## Code

- [code/verify_all.py](code/verify_all.py): the single verifier. Sections `A_*` (symbolic identities), `B_*` (explicit examples, Legendre integrals, cubic sharpness family), `C_*` (all-shift quartic and the Corollary 5.2 identities), `D_*` (persistence for the fixed measure, `--kmax`, default 30), `E_*` (tail family), `F_*` (natural families, the positive theorem and the β = 20 Jacobi example) are exact; `numerical_*` sections use mpmath. `real_root_count` is the exact Sturm diagnostic of Section 8.
- [results/verification.json](results/verification.json): its output.
- [negative_controls.py](negative_controls.py): semantic negative controls.

## Replay

Python 3.13 with the pin in [requirements.txt](requirements.txt) (mpmath 1.3.0; the exact sections need only the standard library):

```sh
sha256sum -c SHA256SUMS                       # or: shasum -a 256 -c SHA256SUMS
python -B -O code/verify_all.py OUT           # all certificates, about 100 s
cmp OUT/verification.json results/verification.json
python -B -O negative_controls.py             # semantic negative controls, about 15 s
```

Results of the release replay: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). Environment: [ENVIRONMENT.txt](ENVIRONMENT.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml).

## Useful next research question

Prove Conjecture 6.4: genuine nonreal poles of δ_k^{(0)}(1; ·) at every k ≥ 3 for μ = δ_1 + 50δ_{1/50} + 10⁻³e^{−t}dt. An effective Rouché bound for Theorem 6.3(2), together with a large-k non-cancellation argument, would close it (Section 9, open problem (i)).
