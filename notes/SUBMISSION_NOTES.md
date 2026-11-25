# Submission notes: "Modular shifts of tadpole Nahm sums"

## What changed from the 9 November draft

**Front matter**
- The manuscript now uses `amsart`, with a new abstract, MSC 2020 codes and keywords.
- Author, affiliation, acknowledgements and data-DOI placeholders are marked in [brackets].

**Introduction (rewritten)**
- New background section: Nahm's problem, Zagier's rank-one classification, Vlasenko–Zwegers, the Calinescu–Milas–Penn conjecture, and the tadpole literature.
- The main results are now stated up front as Theorem 1.2 (classification) and Theorem 1.3 (companions).
- One "Relation to previous work" paragraph replaces the repeated "already known / not a claim of priority" disclaimers.
- New subsections on computations and on notation.
- New remark: Conjecture 1.1 at m = 0 gives exactly Zagier's three rank-one triples with A = 1.

**Draft language removed**
- Deleted: "Consolidated research manuscript", "novelty audit … pending", "independent specialist review has not been completed", "checkpoint".

**Literature added (6 references)**
- Calinescu–Milas–Penn (2016)
- Zagier (2007)
- Nahm (2007)
- Kanade (arXiv:1902.05199)
- Cao–Wang (arXiv:2606.13590)
- Li (arXiv:2606.25866), who independently proved identities (1.12)–(1.14) of Shi–Wang [SW26b], i.e. three of the eight rank-four pairs
- Singular is now cited.

**Notation clashes fixed**
- p = m(2m+1) is now ν; p is kept for the prime.
- The functional L(P) is now Λ(P); L = 2m+3 is unchanged.
- The summation index r in Proposition 5.4 is now i.
- D_m in §2 is now D_ℓ, and the index p in H_k is now g.
- K_j(a) in Proposition 4.4 is now Δ_j(a).
- H_j in §7 is now ℋ_{m,j}.
- The moments M are now 𝔐, so they no longer collide with the set 𝓜_r.
- The Bartlett–Warnaar parameter c is now referred to as "their parameter c".
- A sentence at the start of §4.5 says that x denotes Hall–Littlewood variables from that point on.

**Additions**
- Remark 2.1: consistency checks on the cusp-expansion normalization.
- Remark 7.1: exact q-series evidence at unproved indices, clearly marked as not a proof.
- A Data availability section.

**Mathematics unchanged**
- Numbering is now by section, and all cross-references are automatic.
- Every multi-digit number in the original appears unchanged in the new PDF (checked by machine).

## Independent checks I ran (all passed)
- **Rank 3:** C₀(B) reproduces all six C values, P(B) = 0 at all six points, and E₀ ≡ P(B) identically.
- **Rank 4:** C₄ and R₄ give the eight table values and R₄ = 0 at all eight points. G₁ through G₁₀ vanish at all eight points.
- **Saddle data:** the printed saddle points and covariance M (rank 3), and the saddle points in ranks 4 and 5, satisfy their equations. The cubic and quintic roots lie in the stated intervals.
- **Certificate (6.2) for k = 2:** verified by exact Laurent-polynomial expansion.
- **Identity (4.3), exactly through q¹²:** all (k, z) with m ≤ 3.
- **Identity (4.3), exactly through q¹⁵, at unproved indices:**
  - (8,7) and (9,8) in both sectors;
  - (9,7) and (10,9) with z = 0.
- **Section 5 for m ≤ 4:**
  - Proposition 5.2;
  - C_{j,ℓ} ∈ ℤ[q];
  - C_{j,ℓ} = δ_{jℓ};
  - C_{m+1,ℓ} = (−1)^m A_{m+1,ℓ}.

## Still needed from you before submission
1. **Author details:** name, affiliation, email, acknowledgements and the repository DOI (Zenodo or similar). Many journals require archived data for computer-assisted proofs.
2. **References:** Calinescu–Milas–Penn and Nahm are now confirmed against published listings. Still check the page range of Zagier (2007), and whether Calegari–Garoufalidis–Zagier (arXiv:1712.04887, added in round 2) has a journal version to cite.
3. **Ranks 7 and 9 (strongly recommended):** export cofactor certificates mod p for both containments, Ī ⊆ J̄ and J̄ ⊆ Ī. A standard-library verifier, like the ones you have for ranks 3–5, can then check them without Singular or msolve. Referees are likely to ask for this.
4. **Cusp-expansion formula:** now partly addressed. Remark 2.1 reports a direct numerical test at non-modular shifts in ranks 3–5; the script is `independent_checks/cusp.py`. A written derivation from Vlasenko–Zwegers Theorem 2.3 would still strengthen the paper.
5. **Archive:** list file sizes, SHA-256 checksums and run times in README.txt.
6. **Rank 4 sufficiency:** rows 5–8 depend on unrefereed preprints ([SW26b], with [Li26] covering three of the four). Update the citations if either paper is published before you submit.
7. **Even ranks:** §7 now says ranks 2, 6 and 8 are not treated. Add the reason, for example computational cost or a different structure.
8. **Novelty:** the intro says "To our knowledge … are new". Confirm this after a final arXiv search, since the field moved quickly in 2026.
9. **Journal choice:** fits include Research in Number Theory, Ramanujan J., J. Number Theory, Advances in Mathematics (if Section 5 is strengthened) and Experimental Mathematics (for the computational part). Use the journal's class file if it requires one; the source is plain amsart.

## Round 2 additions
- **Cusp-expansion test:** the Section 2 procedure was compared with Nahm sums evaluated to 50 digits at non-modular shifts. Shifts tested:
  - rank 3: (0.3, −0.7, 0.2) and (−0.4, 0.25, 0.6);
  - rank 4: (0.2, −0.3, 0.1, 0.4) and (0.5, 0.1, −0.6, 0.3);
  - rank 5: (0.3, −0.1, 0.2, −0.5, 0.15).

  Fitted coefficients of ε, ε², ε³ match c₁ − γ₀, κ₂, κ₃ to about 7, 5 and 3 significant digits. At the modular points B = 0 and (1, −1, 1/2) in rank 3, and (1, −1, 1, −1, 1/2) in rank 5, the procedure gives C = −7/80, 9/40 and 67/168, matching your tables.
- **Calegari–Garoufalidis–Zagier** is now cited in the introduction, for their proof of one direction of Nahm's conjecture.
- **Cover letter:** `cover_letter.tex` and `cover_letter.pdf`; fill in the bracketed fields.

## Round 3 additions (strengthening)
- **New Theorem 3.1, rank 2, complete:** exactly 5 modular shifts.
  - Necessity uses the finite-field method plus Lemma 2.3.
  - Sufficiency is proved by hand, via Euler's identity and then Göllnitz–Gordon (two rows) or Lebesgue's identity (three rows).
  - The product normalizations reproduce every C exactly.
- **New Theorem 3.5, rank 6 necessity:** at most 9 shifts.
  - The mod-p quotient has dimension 9, the top parts are finite (dimension 103), and there are 9 explicit rational zeros.
  - Two of the nine are known modular [CTW26]; the other seven are open.
- **New Conjecture 7.1, even ranks:** 2m+3 pairs in rank 2m, plus one exceptional pair in rank 4.
  - It matches ranks 2, 4 and 6 exactly, and the cusp conditions in ranks 8 and 10.
  - In the sector ending in 0, the C formula is the same as in the odd conjecture.
- **New Remark 7.2:** the odd conjecture agrees with the cusp conditions in ranks 11 (18 pairs) and 13 (21 pairs). Unordered 0/1 patterns fail, so the ordering is genuinely forced.
- **New §4.9:** a worked rank-5 example tracing one pair through both halves of the paper.
- **Theorem 1.2** now covers ranks 2, 3, 4, 5, 7 and 9, and the abstract has been updated.
- **References added:** Göllnitz (1967), Gordon (1965) and Andrews (1976), all confirmed against published listings.

### Important for these new results
- **The rank 2 and rank 6 necessity proofs rest on my implementation** (folder `exact_pipeline/` in the zip).
  - It reproduces your rank 3, 4 and 5 results exactly: dimensions 6, 8 and 9, with identical obstruction degrees.
  - Before submission, re-run both eliminations in Singular or msolve, as you did for ranks 7 and 9.
  - Add the `exact_pipeline/` folder to your archive; Appendix B refers to it.
- **Rank 6 sufficiency:** seven of the nine pairs are open. They are not simple infinite products (leading coefficients 5, 9, 5, 6, 15, 20, 30). An even-rank Wronskian formula analogous to (4.3) is the natural next target.
- **Unexplained observation:** for unordered patterns in rank 11, κ₂ and κ₃ take identical values in the z = 0 and z = ½ sectors. This may reflect a symmetry worth explaining.

## Round 4 additions (even-rank Wronskian formula)
- **New Conjecture 7.2: an explicit theta-determinant formula for every even-rank sum.**
  - It mirrors your odd-rank (4.3) with two changes: the lattice carries an alternating sign, and the factor Π vᵢ is dropped.
  - Sector (i) uses half-integer tuples, sectors (ii) and (iii) use integer tuples, and sector (iii) has an extra factor of 2.
  - For B = 0 and (0, …, 0, ½) it is the even-rank Chern–Tran–Wakhare formula.
- **New Proposition 7.3:** Conjecture 7.2 implies modularity with exactly the C values of Conjecture 7.1. The forced C values were checked to agree in every rank 2m with m ≤ 39. Together with Theorem 3.5, it would complete the rank-6 classification.
- **New Remark 7.4, exact q-series verification of every index:**
  - ranks 2 and 4 through q²⁰;
  - rank 6 through q⁴⁰;
  - ranks 8, 10 and 12 through q¹⁵;
  - rank 14 through q¹⁰;
  - the (iv) relation through q²⁰ in ranks 2–8.

  Rank 4 serves as a control, since those sums are already proved modular. The scripts are in `exact_pipeline/even_wronskian/` in the zip.

### Suggested next step (research, not editing)
Prove Conjecture 7.2 for the basic even-rank indices by adapting §4. For example, sector (i) with k = 1 and sector (ii) with j = 1, using the Bartlett–Warnaar identities with the alternating character (CTW already use them in even rank). That would settle at least two more rank-6 shifts and give a general even-rank analogue of Theorem 4.5.

## Round 5: where computation stops
- **Shi–Wang even-rank reduction:** their formula (1.10) was checked numerically in rank 4 and holds with x₀ = 1/x₁. It does not give an even-rank simultaneous shift: a q^{h n₂} term remains. §7.2 now says so, and says that proofs must come from Bartlett–Warnaar at non-principal specializations. The script is `chk5/sw.py`.
- **Chern–Tran–Wakhare** use only the principal specializations, x = 1 with z ∈ {1, q^{1/2}}.
- **Recommendation:** submit now. The classification results (ranks 2, 3, 4, 5, 7, 9, and the rank-6 bound) and the odd-rank companions make a complete paper. The even-rank Wronskian programme (Conjecture 7.2) is a natural follow-up paper and needs genuinely new proofs, not more computation.

## Round 6: even-rank companions proved (new §7.3)
Following Bartlett–Warnaar into even rank:
- **Proposition 7.4: even boundary formulas.** BW Theorem 3.7 holds for any length n. With n = 2m−1 and u = 1, the even-rank sums (0, …, 0, c) become Hall–Littlewood sums in the letters (x̂^±, 1). The b → ∞, c = −w/t limit of BW (5.6b), the type-B identity, gives three cosine-determinant formulas: at t = 1 and t = q with half-integer indices, and at t = w with integer indices. These are the analogues of (4.5)–(4.7).
  - Normalization detail: the frozen letter 1 counts once. It contributes a single (−wt) factor and single cross factors (q x_i^{±1}).
  - Checked numerically to 10⁻¹¹ at generic points in ranks 4 and 6, and against the Nahm sums in ranks 2, 4 and 6.
- **Lemma 7.5: type-B tridiagonal inverse.** Your Lemma 4.7 and Lemma 5.3 go through with one change. Cosine symmetry makes the first diagonal entry 1 + t + t², so the determinant is 1 + t; this explains the factor 1 + w in the t = q formula. For j ≤ m, matrix paths never reach that entry, so the coefficients A_{k,j} are unchanged. Checked numerically to 10⁻¹².
- **Theorem 7.6:** for every m ≥ 2 and 1 ≤ k ≤ min(6, m−1), the even-rank sector (i) Wronskian identity holds.
  - k = 1 uses your contiguity argument.
  - k = 2, …, 6 use your telescoping certificates at u = 1. They are local in the last 2k+1 coordinates, so they apply verbatim in even rank 2m ≥ 2k+1. I verified the resulting reductions exactly in ranks 4–10.
- **Corollary 7.7, rank 6:** four of the nine pairs are now proved modular: (0,…,0), (0,0,0,−1,1,0), (0,−1,1,−1,1,0) and (0,…,0,½).
- **Still open:** sector (i) with k = m, sector (ii) with j ≥ 1, and sectors (iii) and (iv). These need BW-type identities with the unpaired letter u ≠ 1, at u = q^{1/2}. The end of §7.3 says so.

### Please check
- Proposition 7.4(ii) rests on my reconstruction of the b → ∞ limit of BW (5.6b). I read their statement only through summaries, but each consequence matches numerically. Compare the normalization P^B_t with BW's D(x^±; b, c; q) for x = (x̂, 1).
- Theorem 7.6 for k = 3, …, 6 uses your certificate files at u = 1. Confirm that the verifiers don't assume odd rank anywhere. The identities themselves are local, but the file headers may say "rank 2m+1".
- Scripts are in `exact_pipeline/even_wronskian/bw_typeB/` in the zip.

## Round 7: the remaining even sectors
- **x̂ → w doesn't help.** Evaluating the type-B formula at x̂ → w gives no product relation with sectors (i) k = m, (ii) j ≥ 1, (iii) or (iv). I searched quotients against products of (−q;q), (−w;q), (q;q), (q;q²), (−q;q²) and (q²;q²) with exponents up to ±3 and simple constants, at two values of q: nothing. The scripts are in `chk6/`.
- **The obstruction is the unpaired letter w.** The needed alphabet (w^±)^{m−1} ∪ {w} cannot be produced from BW's alphabets (x^±) or (x̂^±, 1). The end of §7.3 now states this precisely. Proving these sectors needs a new Hall–Littlewood identity of BC-type with an unpaired letter. That is a research problem for a follow-up paper, possibly in collaboration with someone from the Warnaar school.
- **New Remark 7.8, even leading coefficients:** C(2m−1, k) − C(2m−1, k−2) in sector (i), and C(2m, j) in sectors (ii)/(iii), the latter including the factor 2 at j = m. Checked for m ≤ 24; the Nahm sums agree in ranks 2–14.
