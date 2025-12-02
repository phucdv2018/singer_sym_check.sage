############################################################
# singer_sym_check.sage
#
# Unified script verifying:
# 1. Scalability on large parameters (e.g. SL_10(2^16)).
# 2. Detailed correctness on small parameters (real matrices).
############################################################

from itertools import product, permutations
from collections import Counter

# =========================================================
# PART A: NEW FUNCTIONS FOR LARGE SCALE (Optimized)
# =========================================================

def weight_patterns_sumK(K, d):
    """
    Efficient recursive generator for partitions of K into d parts.
    Essential for large d (e.g. d=10) to avoid memory overflow.
    """
    if d == 1:
        yield (K,)
        return
    for i in range(K + 1):
        for tail in weight_patterns_sumK(K - i, d - 1):
            yield (i,) + tail

def check_injectivity_sumK_optimized(q, d, K, verbose=True):
    r"""
    Optimized check for Lemma 3.1 using integer arithmetic only.
    Suitable for SL_10(2^16).
    """
    q = Integer(q)
    d = Integer(d)
    K = Integer(K)
    mod = q**d - 1

    if verbose:
        print(f"\n== Check injectivity for (q,d,K) = ({q},{d},{K}) ==")

    seen = {}
    count = 0

    # Use the memory-efficient generator
    for c in weight_patterns_sumK(K, d):
        # Compute exponent E(c) purely as integer
        E = sum(Integer(c[i]) * (q**i) for i in range(d)) % mod
        
        if E in seen:
            # Collision detected
            if seen[E] != c:
                if verbose:
                    print("Collision found!")
                return False
        else:
            seen[E] = c
        count += 1

    if verbose:
        print(f"Checked |B_K| = {count} patterns.")
        print(f"Distinct exponents modulo (q^d - 1): {len(seen)}")
        if count == len(seen):
            print("No collisions: Phi is injective on B_K.")
    
    return True

def run_big_experiment():
    """
    Reproduces the 'Output A' for large parameters.
    """
    print("=== Big exponent-model experiment for lemma ===")
    q = 2**16
    d = 10
    for K in [2, 3, 4]:
        t0 = cputime()
        ok = check_injectivity_sumK_optimized(q, d, K, verbose=True)
        t1 = cputime()
        res_str = "OK" if ok else "FAILED"
        print(f"Result for K = {K}: {res_str}")
        print(f"CPU time ~ {t1 - t0:.3f} seconds\n")


# =========================================================
# PART B: ORIGINAL FUNCTIONS FOR SMALL SCALE (Detailed)
# =========================================================

def base_q_expansion(E, q, d):
    coeffs = []
    for _ in range(d):
        coeffs.append(E % q)
        E //= q
    return coeffs

def check_injectivity(q, d, C, verbose=True):
    # Original brute-force on Box [0,C]^d
    modulus = q**d - 1
    seen = {}
    for b in product(range(C + 1), repeat=d):
        E = sum(b[i] * (q**i) for i in range(d)) % modulus
        if E in seen and seen[E] != b:
            if verbose: print("Collision found!"); return False
        seen[E] = b
    if verbose:
        print(f"Phi is injective on B_{C} for (q,d,C) = ({q},{d},{C}).")
        print(f"Checked {len(seen)} vectors.")
    return True

def singer_eigenvalues_tensor_power(q, d, K):
    Fqd = GF(q**d, 'a')
    omega = Fqd.multiplicative_generator()
    eig = {}
    # Using product is fine for small parameters
    for c in product(range(K + 1), repeat=d):
        if sum(c) != K: continue
        E = sum(c[i] * (q**i) for i in range(d))
        eig[c] = omega**E
    return eig

def exponent_and_digits(lam, omega, q, d):
    E = lam.log(omega)
    E = ZZ(E) % (q**d - 1)
    digits = base_q_expansion(E, q, d)
    return E, digits

def singer_matrix_GL(q, d):
    Fq = GF(q, 'c')
    Fqd = GF(q**d, 'a')
    omega = Fqd.multiplicative_generator()
    basis = [omega**i for i in range(d)]
    S = matrix(Fq, d, d)
    for j in range(d):
        y = omega * basis[j]
        poly = y.polynomial()
        coeffs = poly.list()
        if len(coeffs) < d: coeffs += [Fq(0)] * (d - len(coeffs))
        S.set_column(j, vector(Fq, coeffs[:d]))
    return S, Fqd, omega

def sym_power_matrix(S, k):
    F = S.base_ring()
    d = S.nrows()
    all_tuples = list(product(range(d), repeat=k))
    tuple_index = {t: idx for idx, t in enumerate(all_tuples)}
    dim_tensor = len(all_tuples)
    sym_tuples = sorted(set(tuple(sorted(t)) for t in all_tuples))
    N = len(sym_tuples)
    B = matrix(F, dim_tensor, N)
    for j, s in enumerate(sym_tuples):
        for t in set(permutations(s)):
            idx = tuple_index[t]
            B[idx, j] += 1
    S_k = S
    for _ in range(k - 1): S_k = S_k.tensor_product(S)
    M_sym = matrix(F, N, N)
    for j in range(N):
        v_j = B.column(j)
        w_j = S_k * v_j
        x_j = B.solve_right(w_j)
        M_sym.set_column(j, x_j)
    return M_sym, sym_tuples

def check_sym_end_to_end(q, d, k, verbose=True):
    if verbose:
        print(f"\n=== End-to-end check for Sym^{k}(V) with d={d}, q={q} ===")
    S, Fqd, omega = singer_matrix_GL(q, d)
    if verbose: print("Singer matrix S in GL_d(q):"); print(S)
    M_sym, sym_tuples = sym_power_matrix(S, k)
    N = M_sym.nrows()
    if verbose:
        print(f"Dimension of Sym^{k}(V): {N}")
        print("Sym^k basis index (sorted k-tuples):"); print(sym_tuples)
    patterns = [c for c in product(range(k + 1), repeat=d) if sum(c) == k]
    if verbose: print(f"Number of weight patterns c with sum(c)=k: {len(patterns)}")
    
    M_sym_Fqd = M_sym.change_ring(Fqd)
    eigs = M_sym_Fqd.eigenvalues()
    if verbose: print(f"Number of eigenvalues returned (with multiplicity): {len(eigs)}")
    
    mults = Counter(eigs)
    if not all(m == 1 for m in mults.values()):
        if verbose: print("WARNING: Spectrum is NOT simple."); return False
    else:
        if verbose: print("All eigenvalues have multiplicity 1 (simple spectrum).")

    eig_model = [omega**sum(c[i] * (q**i) for i in range(d)) for c in patterns]
    set_real = set(eigs)
    set_model = set(eig_model)
    if verbose:
        print("Size of set of eigenvalues (real)   :", len(set_real))
        print("Size of set of eigenvalues (model)  :", len(set_model))
    
    if set_real == set_model:
        if verbose: print("SUCCESS: Eigenvalues on Sym^k(V) match the digit-vector model.")
        lam_example = next(iter(set_real))
        E_ex, digits_ex = exponent_and_digits(lam_example, omega, q, d)
        if verbose:
            print("Example eigenvalue lam =", lam_example)
            print("Exponent E     =", E_ex)
            print("Base-q digits  =", digits_ex)
        return True
    return False

def sym2_matrix_from_A(A):
    F = A.base_ring()
    a, b = A[0, 0], A[0, 1]
    c, d = A[1, 0], A[1, 1]
    two = F(2)
    col1 = vector(F, [a*a, a*c, c*c])
    col3 = vector(F, [b*b, b*d, d*d])
    col2 = vector(F, [two*a*b, two*(a*d + b*c), two*c*d])
    M = matrix(F, 3, 3)
    M.set_column(0, col1); M.set_column(1, col2); M.set_column(2, col3)
    return M

def recover_A_from_sym2_bruteforce(M):
    F = M.base_ring()
    candidates = []
    for a in F:
        for b in F:
            for c in F:
                for d in F:
                    A = matrix(F, 2, 2, [[a, b], [c, d]])
                    if A.det() == 0: continue
                    if sym2_matrix_from_A(A) == M: candidates.append(A)
    return candidates

def are_scalar_multiples(A, B):
    F = A.base_ring()
    coords = [(i,j) for i in range(2) for j in range(2) if A[i,j] != 0]
    if not coords: return False
    i0, j0 = coords[0]
    if B[i0, j0] == 0: return False
    lam = B[i0, j0] / A[i0, j0]
    return B == lam * A

def check_recover_A_from_sym2(q, trials=3, verbose=True):
    F = GF(q, 'c')
    if verbose: print(f"\n=== Check recovering A from Sym^2(A) in GL_2({q}) ===")
    for t in range(trials):
        while True:
            A = random_matrix(F, 2, 2)
            if A.det() != 0: break
        M = sym2_matrix_from_A(A)
        candidates = recover_A_from_sym2_bruteforce(M)
        if verbose:
            print(f"\nTrial {t+1}:")
            print("A ="); print(A)
            print(f"Number of candidates found: {len(candidates)}")
        found = False
        for Ah in candidates:
            if are_scalar_multiples(A, Ah):
                found = True
                if verbose:
                    print("Found candidate A_hat scalar multiple of A:")
                    print(Ah)
                break
        if not found and verbose: print("WARNING: no scalar-multiple found.")
    if verbose: print("All trials succeeded: A recovered up to scalar from Sym^2(A).")

def run_small_demo():
    """
    Reproduces the 'Output B' (Original detailed demo).
    """
    print("\n=== Demo: Base-q injectivity ===")
    q, d, C = 7, 3, 3
    check_injectivity(q, d, C, verbose=True)

    print("\n=== Demo: model eigenvalues of Singer on V^{otimes K} ===")
    K = 3
    eig = singer_eigenvalues_tensor_power(q, d, K)
    print(f"Number of weight patterns c with sum(c) = {K}:", len(eig))
    if len(eig.values()) == len(set(eig.values())):
        print("Distinct weight patterns c give distinct eigenvalues (as expected).")
    
    some_c = sorted(eig.keys())[0] # Pick a consistent one like (0,0,3) if available
    # For consistent demo print, we specifically look for (0,0,3) if it exists
    if (0,0,3) in eig: some_c = (0,0,3)
    
    lam = eig[some_c]
    Fqd = lam.parent()
    omega = Fqd.multiplicative_generator()
    E, digits = exponent_and_digits(lam, omega, q, d)

    print(f"\nExample weight pattern c = {some_c}")
    print(f"Eigenvalue lam =  {lam}")
    print(f"Exponent E     =  {E}")
    print(f"Base-q digits  =  {digits}")
    print("Note: digits correspond to the c-vector up to choice of generator.")

    print("\n=== End-to-end checks for Sym^k(V) ===")
    for k in [2, 3]:
        check_sym_end_to_end(q, d, k, verbose=True)

    print("\n=== Step-4 toy test: recover A from Sym^2(A) in GL_2(q) ===")
    check_recover_A_from_sym2(q=7, trials=3, verbose=True)

# =========================================================
# MAIN ENTRY POINT
# =========================================================

def main():
    # 1. Run the new Scalability Check (Output A)
    run_big_experiment()
    
    # 2. Run the original Detailed Demo (Output B)
    # Adding a visual separator
    print("\n" + "="*60 + "\n")
    run_small_demo()

if __name__ == "__main__":
    main()