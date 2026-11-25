# tadpole-nahm-modularity

Data and programs for the paper **Modular shifts of tadpole Nahm sums: classifications, companion identities, and boundary calculus**.

For the tadpole matrix $T_r$, the paper studies the rational shifts $(B,C)$ for which the Nahm sum $F_{r,B,C}$ is a modular function. It contains:
- complete classifications in ranks 2, 3, 4, 5, 7 and 9;
- an upper bound of nine shifts in rank 6;
- companion Wronskian identities in odd rank, and in one even-rank sector;
- a boundary calculus for the Hall–Littlewood side;
- conjectures for all odd and even ranks.

## Layout

| Folder | Contents | Used in a proof? |
|---|---|---|
| `paper/` | LaTeX sources (`main.tex` and section files), the compiled PDF and the cover letter | n/a |
| `computations/elimination/` | Exact cusp polynomials and finite-field eliminations behind Theorem 3.1 (rank 2) and Theorem 3.5 (rank 6) | **yes** |
| `computations/author_data/` | Certificates and data for ranks 3, 4, 5, 7 and 9 and for the companions k = 3,…,6 (see the README there) | **yes** |
| `supplementary/` | Independent checks and numerical evidence. Not cited in the paper and not needed for any proof | no |
| `notes/` | Working notes from the revision | no |

Appendix B of the paper lists exactly the computations that the proofs depend on.

## Reproducing the rank-2 and rank-6 eliminations

These need only Python 3 (standard library) and NumPy for the zero search.

```bash
cd computations/elimination
python3 findpts.py 2 3 2 3   # rank 2: quotient dimension mod p, top parts finite, the five rational zeros
python3 classify.py 3 3      # rank 3: reproduces the known bound of six
python3 run_rank6.py 3       # rank 6: obstruction polynomials, Groebner bases mod p (about 10 minutes); writes r6_N3.pkl
python3 find_rank6_zeros.py 2 3   # rank 6: the nine rational zeros, confirmed exactly
```

`exactcusp.py` carries out the procedure of Section 2.1 in exact arithmetic over Q(t), with t = 1 − cos(π/(r+2)). `gbmodp.py` is a Buchberger implementation over F_p (grevlex), used for both computations required by Lemma 2.3.

## Supplementary material

| Folder | What it checks |
|---|---|
| `supplementary/odd_rank_checks/` | Table data in ranks 3 and 4, the printed k = 2 certificate (6.2), identity (4.3) as exact q-series (including unproved indices), the boundary-coefficient formulas of Section 5, saddle points, and the rank-2 product formulas |
| `supplementary/cusp_numerics/` | Direct numerical test of the Section 2.1 cusp expansion against Nahm sums; cusp-condition tests of the odd and even conjectures in ranks 8–13 |
| `supplementary/even_wronskian/` | Exact q-series tests of the even-rank Wronskian identities (Conjecture 7.2) and the (iv) relation |
| `supplementary/type_b_boundary/` | Numerical confirmation of the type-B boundary formulas (Proposition 7.4), the type-B tridiagonal inverse (Lemma 7.5) and the even-rank Nahm-side reductions; check of the Shi–Wang even-rank reduction; exploration of the remaining even sectors |

## Requirements
- Python ≥ 3.10 and NumPy.
- No computer algebra system is needed for the files in this repository. The rank 7 and 9 computations in `computations/author_data/` use Singular and msolve.
