# Corrections

Errata to release 1.0.0-candidate (DOI 10.5281/zenodo.23586255). The archived PDF, the
`paper/paper.tex` source on this branch (kept byte-identical to the tagged release and
`SHA256SUMS`), the Zenodo deposit and the tag `v1.0.0-candidate` are unchanged.

## 2026-10-23 — §4.1, illustrative remark on internal ratios

The manuscript (§4.1, after the tail-measure family; `paper/paper.tex` line 386–388) says:
"Other internal ratios often fail too (at $k=8$, $b\to0$, five of the seven are below $4$)".

Correct count: **six** of the seven. With $h_j=\binom{8}{j}(b+j)_7$, exact rational arithmetic
gives the internal ratios $h_j^2/(h_{j-1}h_{j+1})$ for $j=2,\dots,7$ at $b=0$ as
28/9, 54/25, 125/66, 11/6, 126/65 and 52/21 (all between 1.83 and 3.12), while the $j=1$
ratio is unbounded as $b\to0$ (about $2\times10^6$ at $b=10^{-6}$). The count is also six of
seven at $b=10^{-6}$ and $b=0.1$.

The remark is not used in any proof; the corrected statement is slightly stronger, so the
robustness point it illustrates stands. No theorem, certificate or evidence state is affected.
Also recorded at <https://evidencepress.org/releases/weniger-delta-nonreal-poles/>.
