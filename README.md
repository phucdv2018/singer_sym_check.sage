# singer_sym_check.sage

SageMath code accompanying the paper:

> A uniform rewriting algorithm for twisted tensor representations of finite classical groups

This script experiments with:
- Base-\(q\) injectivity \(\Phi: B_C \to \mathbb{Z}/(q^d - 1)\mathbb{Z}\),
- Model eigenvalues of a Singer cycle on tensor powers \(V^{\otimes K}\),
- Construction of a real Singer matrix in \(\mathrm{GL}_d(q)\),
- The induced action on \(\mathrm{Sym}^k(V)\) and its eigenvalues,
- End-to-end comparison with the digit-vector model.

## Requirements

- SageMath (version ≥ 9.x recommended).

## Usage

Run the demo from the command line:

```bash
sage singer_sym_check.sage
