# Author data (to be added)

The proofs for ranks 3, 4, 5, 7 and 9 and for the companion identities k = 3,…,6 rely on certificate files from the original computations. They are not in this repository yet. Add them here, keeping the file names quoted in the paper:

- `uniform_k3_certificate.json`, `uniform_k4_certificate.json`
- `uniform_k5_zero_certificate.json`, `uniform_k5_half_certificate.json`
- `uniform_k6_zero_certificate.json`, `uniform_k6_half_certificate.json` (each with field `target_multiplier`)

Also add the certificates and Gröbner data for:
- Lemma A.1 (rank 3): the rational multipliers and the quartics E3, E4;
- Lemma A.2 (rank 4): the obstruction polynomials and elimination certificate;
- Lemma A.3 (rank 5): the obstruction polynomials and elimination certificate;
- Theorems 3.6 and 3.8 (ranks 7 and 9): the obstruction polynomials, the finite-field computations (Singular, msolve) and, for rank 9, the exact saddle data.

Finally, add the programs that check each certificate by exact expansion.
