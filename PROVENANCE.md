# Provenance

* Scholarly creator: **Anonymous**, the standing Evidence Press attribution. No human research contribution is claimed or inferred. There is no ORCID and no affiliation.
* Publisher: Evidence Press (https://evidencepress.org; release page https://evidencepress.org/releases/weniger-delta-nonreal-poles/). Ian Pitchford maintains the repository and carries out publication operations under the repository owner's account. That role is not authorship.
* Production: AI-assisted research. AI research agents produced the manuscript, proofs and computations. Errors were sought in five ways: a second agent checked the proofs and replayed the certificates of the first draft in the same research session (19 October 2026); an internal multi-perspective AI-assisted review of that draft (editor, three specialist perspectives and a devil's advocate, 23 October 2026); an AI-assisted adversarial audit of the same draft by a different model (configured as OpenAI gpt-6-astra, 23 October 2026); a verification re-review of the revised version (23 October 2026); and a second adversarial audit of the revised version by the same external model (23 October 2026). Their substantive points and how each was resolved are in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). The original reports are kept privately.
* Target: R. Borghi, *Is Weniger's transformation capable of simulating the Stieltjes function branch cut?*, Mathematics 14 (2026), no. 2, article 376, doi:10.3390/math14020376 (published 22 January 2026); author preprint https://www.preprints.org/manuscript/202512.2647 (30 December 2025). Section 1 of the manuscript quotes the question of §3.1 and sentences of §3.3 of the published article, with attribution.
* History:
  * A 15-page first draft, dated 19 October 2026, with two verifiers and a bounded novelty search run the same day.
  * The draft was reviewed and audited on 23 October 2026. No reviewer found a mathematical error; the requests concerned positioning, literature, a fixed-measure persistence result, a positive theorem, explicit certificates at k = 5–8, a practical section and exposition.
  * It was revised the same day: retitled, with the persistence theorem (Section 6), the positive theorem (Section 7), Table 3, the n = 1 case of Proposition 5.1, the practical section (Section 8), Tables 1–2, Figure 1 and 13 added references. The two verifiers were merged into `code/verify_all.py`, which also checks the discriminant identity (12).
  * The revised package was replayed from a clean copy and re-reviewed. The re-review found no wrong theorem and listed seven sentence-level corrections. The second adversarial audit, run in parallel on the revised version, found no counterexample to any theorem and raised one claim-scope point and three minor points.
  * A final editor pass made all of these corrections, restored the first draft's all-shift compactness proof as Corollary 5.2, added Remark 7.4 with an exact β = 20 Jacobi certificate (both prompted by the audit), credited Borghi's earlier positive cases, and added the corresponding verifier checks. The final package was replayed byte for byte (15 pages; still 15 with the release front matter).
  * Only the final version is published.
* Publication edits:
  * the reserved DOI and repository URL inserted in the manuscript front matter and in Appendix B, and the PDF rebuilt;
  * README.md replaced by the release README;
  * the release records listed in README.md added, together with `negative_controls.py`, `requirements.txt` and `.github/workflows/replay.yml`, and the licence files `LICENSE` and `LICENSE-CC0`;
  * SHA256SUMS generated.
  No computation script or results file was changed by publication.
