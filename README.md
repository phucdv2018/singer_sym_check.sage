# singer_sym_check.sage

SageMath code accompanying the paper

> *A uniform rewriting algorithm for twisted tensor representations of finite classical groups*

This script provides small-scale computational experiments that support the theoretical results in the paper. It includes:

- Exhaustive checks of the base-\(q\) injectivity map
  \(\Phi: \mathcal{B}_C \to \mathbb{Z}/(q^d - 1)\mathbb{Z}\),
- A model for eigenvalues of a Singer cycle on tensor powers \(V^{\otimes K}\) via digit vectors,
- Construction of an explicit Singer matrix in \(\mathrm{GL}_d(q)\) from multiplication by a generator of \(\mathbb{F}_{q^d}^\times\),
- The induced action of a Singer cycle on \(\mathrm{Sym}^k(V)\) and computation of its eigenvalues over \(\mathbb{F}_{q^d}\),
- End-to-end comparison between the theoretical digit–vector model and the actual spectrum on \(\mathrm{Sym}^k(V)\),
- A toy reconstruction experiment for the rewriting step: recovering a matrix \(A \in \mathrm{GL}_2(q)\) (up to scalar) from its symmetric square \(\mathrm{Sym}^2(A)\).

The code is intended as a collection of sanity checks for small parameters and is **not** optimised for large-scale computations.

## Requirements

- SageMath (version ≥ 9.x recommended).

## Usage

Run the built-in demo from the command line:

```bash
sage singer_sym_check.sage
