# Replay receipt

## Commands

From the root of a fresh extraction, using Python 3.13 with the pin in [requirements.txt](requirements.txt):

```sh
sha256sum -c SHA256SUMS                                   # macOS: shasum -a 256 -c SHA256SUMS
python -B -O code/verify_all.py OUT                       # all certificates (default --kmax 30)
cmp OUT/verification.json results/verification.json
python -B -O negative_controls.py                         # semantic negative controls
```

`verify_all.py` writes only `OUT/verification.json`. `negative_controls.py` works only in temporary directories. Do not pass `results` as `OUT` unless you mean to overwrite the shipped file.

## Runs before publication

| Run | Command | Result |
|---|---|---|
| First-draft verifiers, 19 and 23 October 2026 | `verify.py` and `verify_all_shifts.py` (superseded) | Replayed byte for byte by the producer, by the first-round reviewers and by the first adversarial audit. Superseded by `code/verify_all.py` at revision |
| Revised package (before the final editor pass), 23 October 2026 | `python -B -O code/verify_all.py OUT` | Byte-identical replays by the producer, from a clean copy, by the re-review and by the second adversarial audit (output SHA-256 `2b605165…3cb83f08b`; superseded) |
| Final package, 23 October 2026, two fresh replays by the publication coordinator | `python -B -O code/verify_all.py OUT` | "all exact checks passed"; both byte-identical to `results/verification.json` (SHA-256 `c0cb9b0c1ef6d0870510e55bbd86233b8b53637abc0717f150223366cb3bab84`) |
| Release preparation, final package, 23 October 2026, CPython 3.13.5, mpmath 1.3.0, macOS arm64, from a fresh copy | `python -B -O code/verify_all.py OUT` (about 98 s) | exit 0; byte-identical to `results/verification.json` (SHA-256 `c0cb9b0c…3cb3bab84`) |
| Negative controls, release preparation, final package | `python -B -O negative_controls.py` (about 15 s) | ALL NEGATIVE CONTROLS BEHAVED: baseline passed; 13 corruptions rejected, each for the intended reason; tampered file detected |

## Release-payload replay

This file cannot record a hash of the archive that contains it. The `code/` and `results/` files of the release tree are byte-identical to those of the final package replayed above (see `SHA256SUMS`); the release adds only these records, the licence files and the DOI-bearing manuscript. The frozen archive is replayed again at publication from a fresh extraction (`-B -O`, byte-identical result file, negative controls); that run is part of the publisher's gate record. The repository's CI workflow repeats the SHA-256 check, the full replay and the negative controls on GitHub for the tagged commit.

## What this establishes

Byte-identical replay shows that the shipped output is exactly what the shipped program produces on the stated environment. It does not show that the program implements the manuscript's definitions correctly. That rests on the check of the rational form against the finite-difference definition at 81 points, on the negative controls, and on the independent exact recomputations in the internal reviews (both adversarial audits used their own SymPy code). None of these runs is an independent reproduction by an unaffiliated party. The `numerical_*` sections are floating point; their byte identity is expected only with mpmath 1.3.0.
