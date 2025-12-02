# singer_sym_check.sage

SageMath code accompanying the paper

> *A uniform rewriting algorithm for twisted tensor representations of finite classical groups*

This script provides computational experiments that support the theoretical results in the paper. It ranges from detailed matrix constructions for small parameters to optimized scalability checks for large groups. Key features include:

- **Scalability checks:** Optimized, integer-only verification of the base-$q$ injectivity lemma for large parameters (e.g., $q=2^{16}, d=10$), avoiding explicit field arithmetic.
- **Small-scale exhaustive checks:** Brute-force verification of the map $\Phi: \mathcal{B}_C \to \mathbb{Z}/(q^d - 1)\mathbb{Z}$ for small $d$.
- **Spectral Model:** A model for eigenvalues of a Singer cycle on tensor powers $V^{\otimes K}$ via digit vectors.
- **Real Matrix Construction:** Explicit construction of Singer matrices in $\mathrm{GL}_d(q)$ and their induced action on $\mathrm{Sym}^k(V)$.
- **End-to-end verification:** Comparison between the theoretical digit–vector model and the actual computed spectrum on $\mathrm{Sym}^k(V)$.
- **Toy reconstruction:** A demonstration of the rewriting step: recovering a matrix $A \in \mathrm{GL}_2(q)$ (up to scalar) from its symmetric square $\mathrm{Sym}^2(A)$.

**Note:** While the matrix reconstruction routines are designed for small examples to illustrate correctness, the eigenvalue labeling and injectivity checks are implemented efficiently and verify the theory for large parameters.

## Requirements

- SageMath (version ≥ 9.x recommended).

## Usage

Run the built-in demo from the command line:

```bash
sage singer_sym_check.sage
