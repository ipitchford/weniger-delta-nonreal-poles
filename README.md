# Nonreal poles of Weniger's delta transformation for Stieltjes series: a finite-order classification

**Release 1.0.0-candidate (23 October 2026). Unrefereed research candidate.**
DOI: [10.5281/zenodo.23586255](https://doi.org/10.5281/zenodo.23586255) (this version) ·
concept DOI [10.5281/zenodo.23586254](https://doi.org/10.5281/zenodo.23586254) (all versions) ·
repository <https://github.com/ipitchford/weniger-delta-nonreal-poles> (tag `v1.0.0-candidate`) ·
Evidence Press page <https://evidencepress.org/releases/weniger-delta-nonreal-poles/>.

Checked by internal AI-assisted review only (a multi-perspective review of the first draft, AI-assisted adversarial
audits by a different model of the first draft and of the revised version, and a verification re-review): not peer
review, not independent reproduction, not formal verification. See [STATUS.md](STATUS.md) and [ASSURANCE.md](ASSURANCE.md).

Author: Anonymous. Text and data CC0-1.0; code MIT (see [LICENSES.md](LICENSES.md)).

This package answers a question of R. Borghi, *Is Weniger's transformation capable of simulating the Stieltjes
function branch cut?*, Mathematics 14 (2026), no. 2, 376, doi:[10.3390/math14020376](https://doi.org/10.3390/math14020376):
are all zeros of the denominator of Weniger's delta transformation confined to the negative real axis for every
Stieltjes series? The manuscript is `paper/paper.pdf` (15 pages); its source is `paper/paper.tex`.

## Main results (status in brackets)

Write b = n + β (starting index n ≥ 0, shift β > 0). "Universally safe" means the delta denominator Q_k has only
negative real zeros for every positive measure on [0, ∞) with infinite support and finite moments.

| Order k | Universally safe | Unsafe (some measure gives a genuine nonreal conjugate pair of poles) |
|---|---|---|
| 1, 2 | every n ≥ 0, β > 0 | never |
| 3 | exactly b ≤ b_c = (√33 − 5)/2 ≈ 0.3723 | b > b_c, e.g. every n ≥ 1, and β = 1 |
| ≥ 4 | never | every n ≥ 0, β > 0 |

- The classification above, including genuineness of the nonreal poles [proved; Theorem 4.1]. The witness measure
  depends on (k, n, β).
- Fixed weights against arbitrary Stieltjes moments: Σ h_j z^j/μ_{m+j} is real-rooted for every measure iff
  h_j² ≥ 4 h_{j−1} h_{j+1} for every internal j (Hutchinson–Kurtz sufficiency plus a realisation lemma) [proved; Theorem 3.2].
- Explicit examples at k = 3 and k = 4, and the tail family at k = 5–8 (Table 3) [exact certificates].
- One measure, every shift, order four: for every fixed n some measure gives genuine nonreal poles for every β > 0
  [proved; Corollary 5.2]; the explicit measure λ = δ_1 + 10⁵ δ_{1/100} + 10⁻³ e^{−t} dt does so at n = 0 and n = 1
  [exact certificate of a polynomial identity in β; Proposition 5.1] but not at n = 2.
- One fixed measure μ = δ_1 + 50 δ_{1/50} + 10⁻³ e^{−t} dt, n = 0, β = 1: Q_k(w/k²) converges to an entire function E_μ
  [proved; Theorem 6.1]; nonreal denominator zeros at every sufficiently large order [proved, with one exact
  certificate; Theorem 6.3(2)]; genuine nonreal poles at every order 3 ≤ k ≤ 30 [exact certificate; Theorem 6.3(1)];
  genuine nonreal poles at every k ≥ 3 [open; Conjecture 6.4].
- Positive theorem: real-rootedness at every order when Borghi's P_k is real-rooted, when the reciprocal moments form a
  multiplier sequence, or when μ_{r−1}μ_{r+1} ≥ 4μ_r² [proved; Theorem 7.2]. Covers Γ(ℓ + a) (including the Euler series),
  Γ(νℓ + c) with integer ν, Beta weights with integer exponent, and log-normal moments with σ² ≥ 2 log 2
  [proved; Corollary 7.3]. Borghi proved the superfactorial, Laguerre and integer-difference Jacobi cases at n = 0,
  β = 1 (his §§4.2, 4.3, 4.6); the corollary extends them to every n, β and k with one proof.
- Q_k is real-rooted for every shift iff Borghi's P_k is real-rooted [proved; Remark 7.4]; a non-integer Jacobi weight
  at β = 20 has genuine nonreal poles [exact certificate]. Natural families tested on a finite grid with β ≤ 5
  (975 Sturm counts) [exact certificate on the grid only].
- Pole locations, residues and accuracy (Section 8) [numerical, illustrative].

## Reproduce

Requires Python ≥ 3.9. The exact sections use only the standard library (integers and `fractions.Fraction`). The
`numerical_*` sections use `mpmath` (tested: mpmath 1.3.0 under CPython 3.13.5); without mpmath they are skipped and
the output differs from `results/verification.json`.

```bash
python3 -B -O code/verify_all.py OUTDIR            # all certificates, about 100 s on the release host
cmp OUTDIR/verification.json results/verification.json
python3 -B -O negative_controls.py                 # semantic negative controls, about 15 s
cd paper && pdflatex paper.tex && pdflatex paper.tex && pdflatex paper.tex
```

`verify_all.py` writes `OUTDIR/verification.json` with sorted keys, no timestamps or paths, and floats rounded to
eight significant figures. Every exact check is an explicit `RuntimeError`, so `-O` disables nothing. `--kmax K`
lowers the highest order of the persistence certificate (default 30); the output then differs from the shipped file.

## Files

- `code/verify_all.py`: the single verifier (sections A–F exact, `numerical_*` floating point).
- `results/verification.json`: its output.
- `paper/paper.tex`, `paper/paper.pdf`, `paper/figures/cubic_nonreal_poles.{pdf,png}` (Figure 1; floating-point coordinates).

## Limits

Everything is at fixed finite order. Nothing here says that the delta approximants fail to converge, and the known
positive results for special families stand. For the fixed measure, genuine nonreal poles are proved only for
3 ≤ k ≤ 30; beyond that only nonreal denominator zeros are proved, and k_0 in Theorem 6.3(2) is not effective. The
Weniger–Čížek K_ν family is tested, not proved. No proof assistant was used.

## Release records

- [STATUS.md](STATUS.md): what is proved, certified or numerical, and what is open.
- [ASSURANCE.md](ASSURANCE.md): what has and has not been checked.
- [CLAIMS.json](CLAIMS.json): machine-readable claims, with status and paper location.
- [AI_INDEX.md](AI_INDEX.md): a navigation page for agents.
- [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md): substantive objections from the internal AI-assisted reviews, and their resolution.
- [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md): replay commands and results. [ENVIRONMENT.txt](ENVIRONMENT.txt) and [requirements.txt](requirements.txt): the environment and pinned packages.
- [negative_controls.py](negative_controls.py): semantic negative controls. Run `python -B -O negative_controls.py`; it should reject corrupted claims, code and outputs.
- [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [RESEARCH_GATES.json](RESEARCH_GATES.json): sources, citation checks and the bounded novelty search.
- [PROVENANCE.md](PROVENANCE.md), [LICENSES.md](LICENSES.md), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md): provenance and the licence map.
- `SHA256SUMS`: hashes of every file. Check with `sha256sum -c SHA256SUMS`.

## Citation

Anonymous (2026). *Nonreal poles of Weniger's delta transformation for Stieltjes series: a finite-order
classification* (Version 1.0.0-candidate) [Unrefereed research candidate]. Evidence Press.
https://doi.org/10.5281/zenodo.23586255

Machine-readable metadata: [CITATION.cff](CITATION.cff) and [.zenodo.json](.zenodo.json). Scholarly creator: Anonymous.
Ian Pitchford maintains the repository and publishes it; that role is not authorship.
