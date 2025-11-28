############################################################
# Sage code to experiment with:
#  - Base-q injectivity
#  - Model eigenvalues on tensor powers V^{⊗ K}
#  - Real Singer matrix in GL_d(q)
#  - Real Sym^k(V) action and eigenvalues
#  - Recovering A from Sym^2(A) in a toy example (Step 4 demo).
############################################################

from itertools import product, permutations
from collections import Counter

############################################################
# 1. Base-q expansion
############################################################

def base_q_expansion(E, q, d):
    r"""
    Return the base-q expansion of integer E as a list of length d:

        E = sum_{i=0}^{d-1} c[i] * q^i,  with  0 <= c[i] < q.

    INPUT:
        - E : non-negative integer (0 <= E < q^d is typical)
        - q : prime power (integer >= 2)
        - d : positive integer (dimension)

    OUTPUT:
        - list [c[0], ..., c[d-1]] of digits in base q
    """
    coeffs = []
    for _ in range(d):
        coeffs.append(E % q)
        E //= q
    return coeffs  # c[0], ..., c[d-1]


############################################################
# 2. Check injectivity of Phi : B_C -> Z/(q^d - 1)Z
############################################################

def check_injectivity(q, d, C, verbose=True):
    r"""
    Check injectivity of the map

        Phi : B_C -> Z/(q^d - 1)Z,
        Phi(b_1,...,b_d) = sum_{i=0}^{d-1} b_{i+1} q^i  mod (q^d - 1),

    where
        B_C = { (b_1,...,b_d) in Z_{\ge 0}^d | 0 <= b_i <= C }.

    This is an exhaustive search, so exponential in d.

    INPUT:
        - q : prime power >= 2
        - d : positive integer
        - C : non-negative integer (we are mainly interested in C < q-1)
        - verbose : if True, print result / collisions

    OUTPUT:
        - True  if Phi is injective on B_C
        - False if a collision is found
    """
    modulus = q**d - 1
    seen = {}

    for b in product(range(C + 1), repeat=d):
        E = sum(b[i] * (q**i) for i in range(d)) % modulus
        if E in seen and seen[E] != b:
            if verbose:
                print("Collision found!")
                print("  b =", b, "and", "c =", seen[E], "both map to", E)
            return False
        seen[E] = b

    if verbose:
        print(f"Phi is injective on B_{C} for (q,d,C) = ({q},{d},{C}).")
        print(f"Checked {len(seen)} vectors.")
    return True


############################################################
# 3. Model eigenvalues of a Singer cycle on V^{\otimes K}
############################################################

def singer_eigenvalues_tensor_power(q, d, K):
    r"""
    Model eigenvalues of a Singer cycle on the tensor power V^{\otimes K},

        V = natural d-dimensional module over GF(q).

    We assume that on V ⊗ GF(q^d), a Singer cycle s acts diagonally by
        s.e_i = ell_i * e_i,
    where ell_i = omega^{q^i} for some generator omega of GF(q^d)^*.

    For each vector c = (c_1,...,c_d) in B_K with sum(c) = K, we define
        E(c) = sum_{i=0}^{d-1} c_{i+1} q^i,
        eigenvalue = omega^{E(c)}.

    This matches the abstract formula used in the paper.

    INPUT:
        - q : prime power
        - d : dimension of V
        - K : total "tensor degree"

    OUTPUT:
        - dict eig:
            keys   = tuples c in B_K with sum(c) = K
            values = eigenvalues in GF(q^d)

    NOTE:
        This function verifies the injectivity of the map
            c -> omega^{E(c)} with E(c) = sum c_i q^i,
        for distinct weight patterns c in V^{\otimes K}.

        It does NOT detect possible weight multiplicities in irreducible
        submodules or quotients: in such modules several linearly
        independent vectors may share the same weight pattern c and
        hence the same eigenvalue.
    """
    # Finite field GF(q^d)
    Fqd = GF(q**d, 'a')
    # choose a multiplicative generator of F_{q^d}^*
    omega = Fqd.multiplicative_generator()

    eig = {}
    for c in product(range(K + 1), repeat=d):
        if sum(c) != K:
            continue
        E = sum(c[i] * (q**i) for i in range(d))
        eig_val = omega**E
        eig[c] = eig_val

    return eig


############################################################
# 4. Recover exponent and base-q digits of an eigenvalue
############################################################

def exponent_and_digits(lam, omega, q, d):
    r"""
    Given an eigenvalue lam = omega^E in GF(q^d)^*, return:

        - exponent E  (0 <= E < q^d - 1),
        - its base-q digits [c_0,...,c_{d-1}].

    INPUT:
        - lam   : element of GF(q^d)^*
        - omega : fixed generator of GF(q^d)^*
        - q     : prime power
        - d     : positive integer

    OUTPUT:
        - (E, digits) where digits = [c_0,...,c_{d-1}]
    """
    # discrete log in multiplicative group of GF(q^d)
    E = lam.log(omega)           # lam = omega^E in GF(q^d)^*
    E = ZZ(E) % (q**d - 1)       # normalise exponent modulo q^d - 1

    digits = base_q_expansion(E, q, d)
    return E, digits


############################################################
# 5. Construct a real Singer cycle in GL_d(q)
############################################################

def singer_matrix_GL(q, d):
    r"""
    Construct a Singer cycle in GL_d(q) as the matrix of
    "multiplication by omega" on GF(q^d) viewed as a d-dimensional
    vector space over GF(q).

    INPUT:
        - q : prime power
        - d : positive integer

    OUTPUT:
        - S     : d x d matrix over GF(q) representing multiplication by omega
        - Fqd   : GF(q^d)
        - omega : generator of Fqd = GF(q^d)^*
    """
    Fq = GF(q, 'c')
    Fqd = GF(q**d, 'a')
    # multiplicative generator of F_{q^d}^*
    omega = Fqd.multiplicative_generator()

    # Basis {1, omega, omega^2, ..., omega^{d-1}} of F_{q^d} over F_q
    basis = [omega**i for i in range(d)]

    # Build matrix S: columns are coordinates of omega * basis[j]
    S = matrix(Fq, d, d)
    for j in range(d):
        y = omega * basis[j]        # element in Fqd
        poly = y.polynomial()       # polynomial in Fq[x], y = poly(omega)
        coeffs = poly.list()        # coefficients [a_0, ..., a_m]
        # pad to length d
        if len(coeffs) < d:
            coeffs += [Fq(0)] * (d - len(coeffs))
        elif len(coeffs) > d:
            # should not happen since degree < d, but just in case
            coeffs = coeffs[:d]
        S.set_column(j, vector(Fq, coeffs))

    return S, Fqd, omega


############################################################
# 6. Matrix of Sym^k(V) for a given S in GL_d(q)
############################################################

def sym_power_matrix(S, k):
    r"""
    Given a d x d matrix S over a field F, construct the matrix of the
    induced action on Sym^k(V), where V = F^d is the natural module.

    This is done via:
      - working in V^{⊗ k} with basis indexed by all k-tuples of {0,...,d-1},
      - constructing the symmetric subspace spanned by symmetrised tensors,
      - restricting S^{⊗ k} to this subspace.

    INPUT:
        - S : d x d matrix over a field F
        - k : positive integer

    OUTPUT:
        - M_sym     : N x N matrix over F representing Sym^k(S)
        - sym_tuples: list of sorted k-tuples labelling the basis of Sym^k(V)
    """
    F = S.base_ring()
    d = S.nrows()
    assert S.ncols() == d

    # All k-tuples indexing basis of V^{⊗ k}
    all_tuples = list(product(range(d), repeat=k))
    tuple_index = {t: idx for idx, t in enumerate(all_tuples)}
    dim_tensor = len(all_tuples)

    # Basis of Sym^k(V): sorted k-tuples (nondecreasing)
    sym_tuples = sorted(set(tuple(sorted(t)) for t in all_tuples))
    N = len(sym_tuples)

    # Build matrix B whose columns are symmetrised basis vectors:
    # v_s = sum_{t with sorted(t)=s} e_t  in V^{⊗ k}.
    B = matrix(F, dim_tensor, N)
    for j, s in enumerate(sym_tuples):
        # all distinct permutations of s
        for t in set(permutations(s)):
            idx = tuple_index[t]
            B[idx, j] += 1  # coefficient 1 on each distinct permuted tensor

    # Build S^{⊗ k} as Kronecker product
    S_k = S
    for _ in range(k - 1):
        S_k = S_k.tensor_product(S)

    # Restriction of S_k to Sym^k(V):
    # we want M_sym such that S_k * B = B * M_sym,
    # i.e. for each column j:
    #   w_j = S_k * v_j = B * (col_j of M_sym)
    # solve B * x = w_j for x.
    M_sym = matrix(F, N, N)
    for j in range(N):
        v_j = B.column(j)
        w_j = S_k * v_j
        x_j = B.solve_right(w_j)  # vector in F^N
        M_sym.set_column(j, x_j)

    return M_sym, sym_tuples


############################################################
# 7. End-to-end check for Sym^k(V)
############################################################

def check_sym_end_to_end(q, d, k, verbose=True):
    r"""
    End-to-end test:

      - build a real Singer matrix S in GL_d(q),
      - construct the matrix of S on Sym^k(V),
      - compute eigenvalues on Sym^k(V) over GF(q^d),
      - compare them with the theoretical eigenvalues
        omega^{sum_i c_i q^i} where c ranges over B_k
        with sum(c_i) = k,
      - check that all eigenvalues are simple (multiplicity 1).

    INPUT:
        - q : prime power
        - d : dimension of V
        - k : symmetric power
        - verbose : whether to print details

    OUTPUT:
        - True if the sets of eigenvalues coincide as expected and
          all eigenvalues are simple.
        - False otherwise.
    """
    if verbose:
        print(f"\n=== End-to-end check for Sym^{k}(V) with d={d}, q={q} ===")

    # 1) Real Singer matrix and field
    S, Fqd, omega = singer_matrix_GL(q, d)
    if verbose:
        print("Singer matrix S in GL_d(q):")
        print(S)

    # 2) Matrix on Sym^k(V) over GF(q)
    M_sym, sym_tuples = sym_power_matrix(S, k)
    N = M_sym.nrows()
    if verbose:
        print(f"Dimension of Sym^{k}(V): {N}")
        print("Sym^k basis index (sorted k-tuples):")
        print(sym_tuples)

    # Check dimension equals number of patterns c with sum(c)=k
    patterns = [c for c in product(range(k + 1), repeat=d) if sum(c) == k]
    if verbose:
        print(f"Number of weight patterns c with sum(c)=k: {len(patterns)}")
    if len(patterns) != N and verbose:
        print("WARNING: dim Sym^k(V) != number of weight patterns (unexpected).")

    # 3) Extend scalars to GF(q^d) and compute eigenvalues (with multiplicities)
    M_sym_Fqd = M_sym.change_ring(Fqd)
    eigs = M_sym_Fqd.eigenvalues()
    if verbose:
        print(f"Number of eigenvalues returned (with multiplicity): {len(eigs)}")

    # Check multiplicities explicitly
    mults = Counter(eigs)
    if not all(m == 1 for m in mults.values()):
        if verbose:
            print("WARNING: some eigenvalues have multiplicity > 1.")
        return False
    else:
        if verbose:
            print("All eigenvalues have multiplicity 1 (simple spectrum).")

    # 4) Theoretical eigenvalues from digit vectors c in B_k
    eig_model = [
        omega**sum(c[i] * (q**i) for i in range(d))
        for c in patterns
    ]

    # 5) Compare sets
    set_real = set(eigs)
    set_model = set(eig_model)

    if verbose:
        print("Size of set of eigenvalues (real)   :", len(set_real))
        print("Size of set of eigenvalues (model)  :", len(set_model))

    if set_real == set_model and len(set_real) == len(patterns):
        if verbose:
            print("SUCCESS: Eigenvalues on Sym^k(V) match the digit-vector model.")
        # Example exponent/digits
        lam_example = next(iter(set_real))
        E_ex, digits_ex = exponent_and_digits(lam_example, omega, q, d)
        if verbose:
            print("Example eigenvalue lam =", lam_example)
            print("Exponent E     =", E_ex)
            print("Base-q digits  =", digits_ex)
        return True
    else:
        if verbose:
            print("WARNING: Eigenvalue sets do NOT match as expected.")
        return False


############################################################
# 8. Toy model for Step 4:
#    Recover A from Sym^2(A) in dimension d = 2
############################################################

def sym2_matrix_from_A(A):
    r"""
    Compute the matrix of Sym^2(A) for A in GL_2(F),
    using explicit polynomial formulas.

    Basis of Sym^2(V) is:
        v1 = e_0 ⊗ e_0,
        v2 = e_0 ⊗ e_1 + e_1 ⊗ e_0,
        v3 = e_1 ⊗ e_1.

    For A = [[a,b],[c,d]], in characteristic != 2 we have:

        A.v1 = a^2 v1 + a c v2 + c^2 v3,
        A.v3 = b^2 v1 + b d v2 + d^2 v3,
        A.v2 = 2 a b v1 + 2(a d + b c) v2 + 2 c d v3.

    We build the 3x3 matrix M with these columns.
    """
    F = A.base_ring()
    if F.characteristic() == 2:
        raise ValueError("sym2_matrix_from_A assumes char(F) != 2.")

    a, b = A[0, 0], A[0, 1]
    c, d = A[1, 0], A[1, 1]
    two = F(2)

    col1 = vector(F, [a*a, a*c, c*c])
    col3 = vector(F, [b*b, b*d, d*d])
    col2 = vector(F, [two*a*b, two*(a*d + b*c), two*c*d])

    M = matrix(F, 3, 3)
    M.set_column(0, col1)
    M.set_column(1, col2)
    M.set_column(2, col3)
    return M


def recover_A_from_sym2_bruteforce(M):
    r"""
    Given M = Sym^2(A) for some A in GL_2(F) over a finite field F
    of odd characteristic, try to recover A up to scalar multiples.

    We do this by brute force over all A in GL_2(F), using the explicit
    formula for Sym^2(A) in dimension 2. This is purely experimental and
    intended as a toy model of Step 4 (solving polynomial equations),
    not as an efficient algorithm.

    INPUT:
        - M : 3x3 matrix over a finite field F

    OUTPUT:
        - list of candidate matrices A such that Sym^2(A) = M
    """
    F = M.base_ring()
    q = F.order()
    if F.characteristic() == 2:
        raise ValueError("recover_A_from_sym2_bruteforce assumes char(F) != 2.")

    candidates = []
    for a in F:
        for b in F:
            for c in F:
                for d in F:
                    A = matrix(F, 2, 2, [[a, b], [c, d]])
                    if A.det() == 0:
                        continue
                    M_test = sym2_matrix_from_A(A)
                    if M_test == M:
                        candidates.append(A)

    return candidates


def are_scalar_multiples(A, B):
    r"""
    Test whether 2x2 matrices A and B over a field are equal up to a nonzero scalar.

    Return True if there exists lambda != 0 in the base field such that B = lambda * A.
    """
    F = A.base_ring()
    assert B.base_ring() == F
    # find first nonzero entry of A
    idx = None
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            if A[i, j] != 0:
                idx = (i, j)
                break
        if idx is not None:
            break
    if idx is None:
        # A is the zero matrix, but we don't consider this case in GL_2
        return False

    i0, j0 = idx
    if B[i0, j0] == 0:
        return False
    lam = B[i0, j0] / A[i0, j0]

    # check all entries
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            if lam * A[i, j] != B[i, j]:
                return False
    return True


def check_recover_A_from_sym2(q, trials=5, verbose=True):
    r"""
    Experimental check of recovering A from Sym^2(A) in GL_2(q).

    For each trial:
        - pick a random A in GL_2(q),
        - compute M = Sym^2(A) via explicit formulas,
        - run recover_A_from_sym2_bruteforce(M),
        - check that among the candidates there is some A_hat that
          is a scalar multiple of A.

    This demonstrates, in a small 2-dimensional case, the idea of
    Step 4: reconstructing the underlying matrix from its action on
    a polynomial module by solving a system of polynomial equations.
    """
    F = GF(q, 'c')
    if F.characteristic() == 2:
        raise ValueError("check_recover_A_from_sym2 is intended for odd q.")

    if verbose:
        print(f"\n=== Check recovering A from Sym^2(A) in GL_2({q}) ===")

    success = True
    for t in range(trials):
        # random invertible A
        while True:
            A = random_matrix(F, 2, 2)
            if A.det() != 0:
                break

        M = sym2_matrix_from_A(A)
        candidates = recover_A_from_sym2_bruteforce(M)
        if verbose:
            print(f"\nTrial {t+1}:")
            print("A =")
            print(A)
            print(f"Number of candidates found: {len(candidates)}")

        found = False
        for Ah in candidates:
            if are_scalar_multiples(A, Ah):
                found = True
                if verbose:
                    print("Found candidate A_hat scalar multiple of A:")
                    print(Ah)
                break

        if not found:
            success = False
            if verbose:
                print("WARNING: no scalar-multiple candidate found for this A.")

    if success and verbose:
        print("All trials succeeded: A recovered up to scalar from Sym^2(A).")
    return success


############################################################
# 9. Simple demo / sanity checks
############################################################

def demo():
    """
    Run a few small demonstrations:

      - Check base-q injectivity for some (q,d,C).
      - Compute model eigenvalues for V^{⊗ K} and verify injectivity.
      - Run end-to-end checks for Sym^k(V) with real Singer matrices.
      - Run a toy Step-4 test: recover A from Sym^2(A) in dimension 2.
    """
    print("=== Demo: Base-q injectivity ===")
    q = 7
    d = 3
    C = 3   # For the paper, we care about C = K with K < q-1.
    check_injectivity(q, d, C, verbose=True)

    print("\n=== Demo: model eigenvalues of Singer on V^{⊗ K} ===")
    K = 3
    eig = singer_eigenvalues_tensor_power(q, d, K)
    print(f"Number of weight patterns c with sum(c) = {K}:", len(eig))

    # Check that distinct c's give distinct eigenvalues in this model:
    vals = list(eig.values())
    if len(vals) == len(set(vals)):
        print("Distinct weight patterns c give distinct eigenvalues (as expected).")
    else:
        print("Some eigenvalues coincide (this should not happen if K < q-1).")

    # Take one eigenvalue and recover its exponent/digits
    some_c = next(iter(eig.keys()))
    lam = eig[some_c]
    Fqd = lam.parent()
    omega = Fqd.multiplicative_generator()
    E, digits = exponent_and_digits(lam, omega, q, d)

    print("\nExample weight pattern c =", some_c)
    print("Eigenvalue lam = ", lam)
    print("Exponent E     = ", E)
    print("Base-q digits  = ", digits)
    print("Note: digits correspond to the c-vector up to choice of generator.")

    # End-to-end tests for Sym^k(V) with real Singer matrices
    print("\n=== End-to-end checks for Sym^k(V) ===")
    for k in [2, 3]:
        ok = check_sym_end_to_end(q, d, k, verbose=True)
        print(f"End-to-end check for Sym^{k}(V):", "OK" if ok else "FAILED")

    # Toy Step-4 test: recover A from Sym^2(A) in GL_2(q)
    print("\n=== Step-4 toy test: recover A from Sym^2(A) in GL_2(q) ===")
    check_recover_A_from_sym2(q=7, trials=3, verbose=True)


# If you run this file as a script in Sage, run the demo.
if __name__ == "__main__":
    demo()
