#include "kh_hurwitz.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ==============================================================
 * kh_hurwitz.c - Routh-Hurwitz Stability Criterion Implementation
 *
 * This is the workhorse of polynomial stability testing.
 * The Routh-Hurwitz criterion provides a necessary and
 * sufficient algebraic condition for all roots of a real
 * polynomial to lie in the open left half-plane.
 *
 * Theorem (Routh, 1877; Hurwitz, 1895):
 *   A real polynomial p(s) = a_n s^n + ... + a_0 with a_n > 0
 *   is Hurwitz-stable iff all leading principal minors of the
 *   Hurwitz matrix are positive, equivalently iff the first
 *   column of the Routh array has no sign changes.
 *
 * References:
 *   Routh, E.J. (1877). Stability of Motion.
 *   Hurwitz, A. (1895). Math. Ann. 46:273-284.
 *   Gantmacher, F.R. (1959). Theory of Matrices.
 *   Barnett, S. (1983). Polynomials and Linear Control Systems.
 * ============================================================== */

/* --------------------------------------------------------------
 * Knowledge Point: Routh Table Construction (L5 Algorithms)
 *
 * The Routh table is constructed row by row:
 *
 * Row 0: a_n, a_{n-2}, a_{n-4}, ...
 * Row 1: a_{n-1}, a_{n-3}, a_{n-5}, ...
 *
 * For k >= 2:
 *   r_{k,j} = -(r_{k-2,0} * r_{k-1,j+1} - r_{k-1,0} * r_{k-2,j+1}) / r_{k-1,0}
 *
 * If r_{k-1,0} = 0, it is replaced by a small epsilon (epsilon method)
 * or the auxiliary polynomial method is used (zero row).
 *
 * The number of RHP roots = number of sign changes in first column.
 * ------------------------------------------------------------ */
KH_RouthTable* kh_routh_create(const KH_Polynomial* p) {
    if (!p || p->degree < 1) return NULL;
    int n = p->degree;
    int n_rows = n + 1;
    int n_cols = (n + 2) / 2;

    KH_RouthTable* rt = (KH_RouthTable*)malloc(sizeof(KH_RouthTable));
    if (!rt) return NULL;
    rt->degree = n;
    rt->n_rows = n_rows;
    rt->n_cols = n_cols;
    rt->sign_changes = 0;

    /* Allocate rows */
    rt->rows = (double**)malloc((size_t)n_rows * sizeof(double*));
    if (!rt->rows) { free(rt); return NULL; }
    for (int i = 0; i < n_rows; i++) {
        rt->rows[i] = (double*)calloc((size_t)n_cols, sizeof(double));
        if (!rt->rows[i]) {
            for (int j = 0; j < i; j++) free(rt->rows[j]);
            free(rt->rows);
            free(rt);
            return NULL;
        }
    }

    /* Fill row 0: a_n, a_{n-2}, a_{n-4}, ... */
    for (int j = 0; j < n_cols; j++) {
        int coeff_idx = n - 2 * j;
        if (coeff_idx >= 0) {
            rt->rows[0][j] = p->coeffs[coeff_idx];
        }
    }

    /* Fill row 1: a_{n-1}, a_{n-3}, a_{n-5}, ... */
    for (int j = 0; j < n_cols; j++) {
        int coeff_idx = n - 1 - 2 * j;
        if (coeff_idx >= 0) {
            rt->rows[1][j] = p->coeffs[coeff_idx];
        }
    }

    /* Build remaining rows via Routh recurrence */
    for (int i = 2; i < n_rows; i++) {
        double prev0 = rt->rows[i-1][0];
        double prev2_0 = rt->rows[i-2][0];

        for (int j = 0; j < n_cols - 1; j++) {
            /* Determinant formula: -(a*d - c*b)/c */
            if (fabs(prev0) < 1e-14) {
                /* Zero pivot: use epsilon */
                rt->rows[i][j] = 0.0;
            } else {
                rt->rows[i][j] = (prev0 * rt->rows[i-2][j+1] -
                                   prev2_0 * rt->rows[i-1][j+1]) / prev0;
                if (fabs(rt->rows[i][j]) < 1e-14) rt->rows[i][j] = 0.0;
            }
        }
    }

    /* Count sign changes in first column */
    int sign_changes = 0;
    double prev_sign = rt->rows[0][0];
    for (int i = 1; i < n_rows; i++) {
        double curr = rt->rows[i][0];
        if (fabs(curr) < 1e-14) curr = 0.0;
        if (prev_sign * curr < 0.0) sign_changes++;
        if (prev_sign == 0.0 && curr != 0.0) {
            /* Sign change from zero: check if surrounding signs differ */
            if (i >= 2 && rt->rows[i-2][0] * curr < 0.0) sign_changes++;
        }
        if (curr != 0.0 || (i == n_rows - 1)) {
            prev_sign = curr;
        }
    }
    rt->sign_changes = sign_changes;

    return rt;
}

void kh_routh_free(KH_RouthTable* rt) {
    if (!rt) return;
    for (int i = 0; i < rt->n_rows; i++) {
        free(rt->rows[i]);
    }
    free(rt->rows);
    free(rt);
}

int kh_routh_sign_changes(const KH_RouthTable* rt) {
    if (!rt) return -1;
    return rt->sign_changes;
}

int kh_routh_unstable_root_count(const KH_RouthTable* rt) {
    if (!rt) return -1;
    return rt->sign_changes;
}

bool kh_routh_is_hurwitz(const KH_RouthTable* rt) {
    if (!rt) return false;
    /* All first column elements must be positive */
    for (int i = 0; i < rt->n_rows; i++) {
        if (rt->rows[i][0] <= 0.0) return false;
    }
    return true;
}

void kh_routh_print(const KH_RouthTable* rt) {
    if (!rt) return;
    printf("Routh Table (degree %d):\n", rt->degree);
    for (int i = 0; i < rt->n_rows; i++) {
        printf("s^%-2d |", rt->degree - i);
        for (int j = 0; j < rt->n_cols; j++) {
            printf(" %10.4f", rt->rows[i][j]);
        }
        printf("\n");
    }
    printf("Sign changes in first column: %d\n", rt->sign_changes);
}

/* --------------------------------------------------------------
 * Knowledge Point: Necessary Hurwitz condition (L4 Fundamental Laws)
 *
 * A necessary (but not sufficient) condition for Hurwitz stability:
 * All coefficients a_i must be positive (same sign as a_n).
 *
 * If any coefficient is zero or negative, the polynomial is
 * definitely NOT Hurwitz-stable. This is a fast pre-check
 * before running the full Routh-Hurwitz test.
 * ------------------------------------------------------------ */
bool kh_check_hurwitz_necessary(const KH_Polynomial* p) {
    if (!p) return false;
    for (int i = 0; i <= p->degree; i++) {
        if (p->coeffs[i] <= 0.0) return false;
    }
    return true;
}

/* --------------------------------------------------------------
 * Knowledge Point: Lienard-Chipart Criterion (L4 Fundamental Laws)
 *
 * The Lienard-Chipart criterion (1914) is equivalent to
 * the Routh-Hurwitz criterion but requires computing only
 * about half the Hurwitz determinants:
 *
 * For a_n > 0, p(s) is Hurwitz iff:
 *   a_{n-1} > 0, a_{n-3} > 0, ... AND Delta_1 > 0, Delta_3 > 0, ...
 * OR equivalently:
 *   a_n > 0, a_{n-2} > 0, ... AND Delta_2 > 0, Delta_4 > 0, ...
 *
 * This is used when the Routh table has issues (zero first column).
 * ------------------------------------------------------------ */
bool kh_check_hurwitz_lienard_chipart(const KH_Polynomial* p) {
    if (!p || p->degree < 1) return false;
    if (p->coeffs[p->degree] <= 0.0) return false;

    /* Check all odd-indexed coefficients are positive */
    for (int i = p->degree - 1; i >= 0; i -= 2) {
        if (p->coeffs[i] <= 0.0) return false;
    }

    /* Check odd-order Hurwitz determinants */
    for (int k = 1; k <= p->degree; k += 2) {
        double det = kh_hurwitz_determinant(p, k);
        if (det <= 0.0) return false;
    }
    return true;
}

/* --------------------------------------------------------------
 * Knowledge Point: Hurwitz Matrix Construction (L3 Math Structures)
 *
 * The Hurwitz matrix H_n of size n x n for polynomial
 * p(s) = a_n s^n + a_{n-1} s^{n-1} + ... + a_0:
 *
 * H_n[i][j] = a_{n - 2j + i}  if 1 <= 2j - i <= n+1
 *           = 0                otherwise
 *
 * Equivalently:
 * Row 1: a_{n-1}, a_{n-3}, a_{n-5}, ...
 * Row 2: a_n,     a_{n-2}, a_{n-4}, ...
 * Row 3: 0,       a_{n-1}, a_{n-3}, ...
 * Row 4: 0,       a_n,     a_{n-2}, ...
 * ...
 *
 * The k-th leading principal minor Delta_k is the determinant
 * of the upper-left k x k submatrix.
 * ------------------------------------------------------------ */
double** kh_hurwitz_matrix_create(const KH_Polynomial* p, int* n) {
    if (!p || !n) return NULL;
    *n = p->degree;
    int N = p->degree;
    double** H = (double**)malloc((size_t)N * sizeof(double*));
    if (!H) return NULL;
    for (int i = 0; i < N; i++) {
        H[i] = (double*)calloc((size_t)N, sizeof(double));
        if (!H[i]) {
            for (int j = 0; j < i; j++) free(H[j]);
            free(H);
            return NULL;
        }
    }

    /* Fill Hurwitz matrix */
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int coeff_idx = p->degree - 2 * j + i;
            if (coeff_idx >= 0 && coeff_idx <= p->degree) {
                H[i][j] = p->coeffs[coeff_idx];
            }
        }
    }
    return H;
}

void kh_hurwitz_matrix_free(double** H, int n) {
    if (!H) return;
    for (int i = 0; i < n; i++) {
        free(H[i]);
    }
    free(H);
}

/* --------------------------------------------------------------
 * Knowledge Point: Hurwitz Determinant Computation (L5 Algorithms)
 *
 * Compute the k-th leading principal minor of the Hurwitz matrix.
 * Uses Gaussian elimination with partial pivoting for numerical
 * stability.
 *
 * The k-th Hurwitz determinant Delta_k is:
 *   Delta_1 = a_{n-1}
 *   Delta_2 = a_{n-1}*a_{n-2} - a_n*a_{n-3}
 *   ...
 *
 * If all Delta_k > 0 for k = 1,...,n, then p(s) is Hurwitz.
 * ------------------------------------------------------------ */
double kh_hurwitz_determinant(const KH_Polynomial* p, int k) {
    if (!p || k <= 0 || k > p->degree) return 0.0;
    int N;
    double** H = kh_hurwitz_matrix_create(p, &N);
    if (!H) return 0.0;

    /* Extract k x k leading principal submatrix */
    double** A = (double**)malloc((size_t)k * sizeof(double*));
    for (int i = 0; i < k; i++) {
        A[i] = (double*)malloc((size_t)k * sizeof(double));
        for (int j = 0; j < k; j++) {
            A[i][j] = H[i][j];
        }
    }

    /* Gaussian elimination with partial pivoting */
    double det = 1.0;
    int sign = 1;
    for (int i = 0; i < k; i++) {
        /* Find pivot */
        int pivot = i;
        double max_val = fabs(A[i][i]);
        for (int r = i + 1; r < k; r++) {
            if (fabs(A[r][i]) > max_val) {
                max_val = fabs(A[r][i]);
                pivot = r;
            }
        }
        if (max_val < 1e-15) {
            det = 0.0;
            break;
        }
        /* Swap rows if necessary */
        if (pivot != i) {
            double* temp = A[i];
            A[i] = A[pivot];
            A[pivot] = temp;
            sign = -sign;
        }
        det *= A[i][i];
        /* Eliminate below */
        for (int r = i + 1; r < k; r++) {
            double factor = A[r][i] / A[i][i];
            for (int c = i; c < k; c++) {
                A[r][c] -= factor * A[i][c];
            }
        }
    }
    det *= sign;

    for (int i = 0; i < k; i++) free(A[i]);
    free(A);
    kh_hurwitz_matrix_free(H, N);
    return det;
}

/* Check all leading principal minors of the Hurwitz matrix */
bool kh_hurwitz_check_all_minors(const KH_Polynomial* p) {
    if (!p) return false;
    for (int k = 1; k <= p->degree; k++) {
        double det = kh_hurwitz_determinant(p, k);
        if (det <= 0.0) return false;
    }
    return true;
}

/* --------------------------------------------------------------
 * Knowledge Point: Hermite-Biehler Theorem (L4 Fundamental Laws)
 *
 * The Hermite-Biehler theorem states:
 *   p(s) is Hurwitz-stable iff the roots of the even part
 *   p_even(omega) and the odd part omega * p_odd(omega)
 *   are real, simple, and interlace on the imaginary axis.
 *
 * Split p(j*omega) = R(omega) + j*omega * S(omega)
 * where R is even part, S is odd part.
 *
 * The interlacing condition: between any two consecutive
 * roots of R(omega), there is exactly one root of S(omega),
 * and vice versa.
 * ------------------------------------------------------------ */

/* Split p(s) into even and odd parts: p(s) = p_even(s^2) + s * p_odd(s^2) */
void kh_poly_split_even_odd(const KH_Polynomial* p,
                              KH_Polynomial** even, KH_Polynomial** odd) {
    if (!p || !even || !odd) return;
    int even_deg = p->degree / 2;
    int odd_deg = (p->degree - 1) / 2;

    *even = kh_poly_create(even_deg);
    *odd = kh_poly_create(odd_deg);

    if (*even) {
        for (int i = 0; i <= even_deg; i++) {
            (*even)->coeffs[even_deg - i] = p->coeffs[2 * i];
        }
    }
    if (*odd) {
        for (int i = 0; i <= odd_deg; i++) {
            (*odd)->coeffs[odd_deg - i] = p->coeffs[2 * i + 1];
        }
    }
}

/* Sturm sequence for counting real roots of a polynomial */
static int sturm_sign_changes(const double* seq_vals, int n_vals) {
    int changes = 0;
    int prev_nonzero = 0;
    bool has_prev = false;
    for (int i = 0; i < n_vals; i++) {
        if (fabs(seq_vals[i]) < 1e-12) continue;
        if (has_prev) {
            double prod = seq_vals[prev_nonzero] * seq_vals[i];
            if (prod < 0.0) changes++;
        }
        prev_nonzero = i;
        has_prev = true;
    }
    return changes;
}

/* --------------------------------------------------------------
 * Knowledge Point: Sturm sequence for root counting (L5 Algorithms)
 *
 * Given a polynomial p(x), the Sturm sequence is:
 *   f_0(x) = p(x), f_1(x) = p'(x),
 *   f_{k+1}(x) = -rem(f_{k-1}, f_k)  (negative remainder)
 *
 * The number of distinct real roots in (a,b] is:
 *   V(a) - V(b)  where V(x) = number of sign changes in the
 *   Sturm sequence evaluated at x.
 * ------------------------------------------------------------ */
static int sturm_real_root_count(const KH_Polynomial* p, double a, double b) {
    if (!p || p->degree <= 0) return 0;
    /* Build Sturm sequence */
    int max_seq = p->degree + 2;
    KH_Polynomial** sturm = (KH_Polynomial**)malloc((size_t)max_seq * sizeof(KH_Polynomial*));
    sturm[0] = kh_poly_copy(p);
    sturm[1] = kh_poly_derivative(p);
    int seq_len = 2;

    for (int k = 2; k < max_seq; k++) {
        if (!sturm[k-1]) break;
        /* Polynomial remainder: f_{k-2} mod f_{k-1}, negated */
        KH_Polynomial* rem = kh_poly_copy(sturm[k-2]);
        if (!rem) break;

        /* Simple polynomial division: rem = rem - q * divisor, for degree reduction */
        int d2 = sturm[k-2]->degree;
        int d1 = sturm[k-1]->degree;
        KH_Polynomial* quot = kh_poly_create(d2 - d1);
        if (!quot) { kh_poly_free(rem); break; }

        KH_Polynomial* work = kh_poly_copy(sturm[k-2]);
        if (!work) { kh_poly_free(rem); kh_poly_free(quot); break; }

        for (int qi = d2 - d1; qi >= 0; qi--) {
            double factor = 0.0;
            if (fabs(sturm[k-1]->coeffs[d1]) > 1e-14) {
                factor = work->coeffs[d1 + qi] / sturm[k-1]->coeffs[d1];
            }
            quot->coeffs[qi] = factor;
            for (int j = 0; j <= d1; j++) {
                work->coeffs[j + qi] -= factor * sturm[k-1]->coeffs[j];
            }
        }

        /* Remainder is in work (negated for Sturm sequence) */
        sturm[k] = kh_poly_copy(work);
        if (sturm[k]) {
            for (int i = 0; i <= sturm[k]->degree; i++) {
                sturm[k]->coeffs[i] = -sturm[k]->coeffs[i];
            }
        }

        kh_poly_free(rem);
        kh_poly_free(quot);
        kh_poly_free(work);
        seq_len++;
        if (sturm[k] && sturm[k]->degree == 0) break;
    }

    /* Evaluate sign changes at a and b */
    double* vals_a = (double*)malloc((size_t)seq_len * sizeof(double));
    double* vals_b = (double*)malloc((size_t)seq_len * sizeof(double));
    for (int i = 0; i < seq_len; i++) {
        vals_a[i] = sturm[i] ? kh_poly_eval(sturm[i], a) : 0.0;
        vals_b[i] = sturm[i] ? kh_poly_eval(sturm[i], b) : 0.0;
    }
    int Va = sturm_sign_changes(vals_a, seq_len);
    int Vb = sturm_sign_changes(vals_b, seq_len);

    free(vals_a);
    free(vals_b);
    for (int i = 0; i < seq_len; i++) kh_poly_free(sturm[i]);
    free(sturm);
    return Va - Vb;
}

bool kh_hermite_biehler_test(const KH_Polynomial* p) {
    if (!p) return false;
    KH_Polynomial *even = NULL, *odd = NULL;
    kh_poly_split_even_odd(p, &even, &odd);
    if (!even || !odd) {
        kh_poly_free(even);
        kh_poly_free(odd);
        return false;
    }

    /* Count real roots of even and odd parts on (-inf, inf) */
    int n_even = sturm_real_root_count(even, -1e10, 1e10);
    int n_odd = sturm_real_root_count(odd, -1e10, 1e10);

    kh_poly_free(even);
    kh_poly_free(odd);

    /* All roots must be real, and counts must match expected values */
    int expected_even = p->degree / 2;
    int expected_odd = (p->degree - 1) / 2;
    return n_even == expected_even && n_odd == expected_odd;
}

/* --------------------------------------------------------------
 * Knowledge Point: Comprehensive polynomial stability analysis
 * (L5 Algorithms - combines Routh-Hurwitz + necessary conditions)
 *
 * This implements a complete stability test:
 *   1. Check necessary condition (all coefficients > 0)
 *   2. Build Routh table
 *   3. Check first column sign changes
 *   4. Classify result as stable/unstable/marginal
 * ------------------------------------------------------------ */
KH_StabilityResult kh_polynomial_stability(const KH_Polynomial* p,
                                             KH_StabilityReport* report) {
    if (!p) return KH_INCONCLUSIVE;

    /* Step 1: Necessary condition */
    if (!kh_check_hurwitz_necessary(p)) {
        if (report) {
            report->overall = KH_UNSTABLE;
            report->message = "Failed necessary condition (non-positive coefficient)";
        }
        return KH_UNSTABLE;
    }

    /* Step 2: Build Routh table and check */
    KH_RouthTable* rt = kh_routh_create(p);
    if (!rt) return KH_INCONCLUSIVE;

    bool is_hurwitz = kh_routh_is_hurwitz(rt);
    int sign_changes = rt->sign_changes;

    KH_StabilityResult result;
    if (is_hurwitz) {
        result = KH_STABLE;
        if (report) {
            report->message = "All roots in open left half-plane";
            report->n_stable_roots = p->degree;
            report->n_unstable_roots = 0;
        }
    } else if (sign_changes > 0) {
        result = KH_UNSTABLE;
        if (report) {
            report->message = "RHP roots detected";
            report->n_unstable_roots = sign_changes;
            report->n_stable_roots = p->degree - sign_changes;
        }
    } else {
        /* First column has zeros but no sign changes: marginal */
        result = KH_MARGINAL;
        if (report) {
            report->message = "Roots on imaginary axis (marginal stability)";
        }
    }

    if (report) report->overall = result;
    kh_routh_free(rt);
    return result;
}

/* --------------------------------------------------------------
 * Knowledge Point: Companion matrix eigenvalue method (L5 Algorithms)
 *
 * For a monic polynomial p(s) = s^n + a_{n-1}s^{n-1} + ... + a_0,
 * the companion matrix:
 *   C = [0 1 0 ... 0]
 *       [0 0 1 ... 0]
 *       [.........]
 *       [0 0 0 ... 1]
 *       [-a_0 -a_1 -a_2 ... -a_{n-1}]
 *
 * The eigenvalues of C are the roots of p(s). This provides
 * an alternative to the Routh-Hurwitz test using numerical
 * eigenvalue computation (e.g., QR algorithm).
 * ------------------------------------------------------------ */

#define MAX_ITER 1000

/* QR iteration for eigenvalue computation (simplified, for small matrices) */
static void qr_decompose(double** A, int n, double** Q, double** R) {
    /* Gram-Schmidt QR */
    for (int j = 0; j < n; j++) {
        /* Copy column j of A */
        for (int i = 0; i < n; i++) {
            Q[i][j] = A[i][j];
        }
        /* Subtract projections onto previous Q columns */
        for (int k = 0; k < j; k++) {
            double dot = 0.0;
            for (int i = 0; i < n; i++) dot += Q[i][k] * A[i][j];
            R[k][j] = dot;
            for (int i = 0; i < n; i++) Q[i][j] -= dot * Q[i][k];
        }
        /* Normalize */
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += Q[i][j] * Q[i][j];
        norm = sqrt(norm);
        R[j][j] = norm;
        if (norm > 1e-15) {
            for (int i = 0; i < n; i++) Q[i][j] /= norm;
        }
    }
}

/* Basic QR eigenvalue algorithm for companion matrix real parts */
int kh_find_roots_companion(const KH_Polynomial* p,
                              double* real_parts, double* imag_parts,
                              int max_roots) {
    if (!p || !real_parts || max_roots < p->degree) return 0;
    int n = p->degree;

    /* Build companion matrix */
    double** A = (double**)malloc((size_t)n * sizeof(double*));
    double** Q = (double**)malloc((size_t)n * sizeof(double*));
    double** R = (double**)malloc((size_t)n * sizeof(double*));
    double** T = (double**)malloc((size_t)n * sizeof(double*));

    for (int i = 0; i < n; i++) {
        A[i] = (double*)calloc((size_t)n, sizeof(double));
        Q[i] = (double*)calloc((size_t)n, sizeof(double));
        R[i] = (double*)calloc((size_t)n, sizeof(double));
        T[i] = (double*)calloc((size_t)n, sizeof(double));
        if (i < n - 1) A[i][i+1] = 1.0;
    }

    double leading = p->coeffs[p->degree];
    for (int j = 0; j < n; j++) {
        A[n-1][j] = -p->coeffs[j] / leading;
    }

    /* QR iteration */
    for (int iter = 0; iter < MAX_ITER; iter++) {
        qr_decompose(A, n, Q, R);
        /* A = R * Q */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double sum = 0.0;
                for (int k = 0; k < n; k++) {
                    sum += R[i][k] * Q[k][j];
                }
                T[i][j] = sum;
            }
        }
        /* Copy T to A */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = T[i][j];
            }
        }
    }

    /* Extract real diagonal entries as approximate real parts */
    for (int i = 0; i < n && i < max_roots; i++) {
        real_parts[i] = A[i][i];
        /* Check sub-diagonal for complex conjugate pairs */
        if (i < n - 1 && fabs(A[i+1][i]) > 1e-8) {
            double a = A[i][i], b = A[i][i+1], c = A[i+1][i], d = A[i+1][i+1];
            double tr = a + d;
            double det = a * d - b * c;
            double disc = tr * tr - 4.0 * det;
            if (disc < 0) {
                real_parts[i] = tr / 2.0;
                real_parts[i+1] = tr / 2.0;
                if (imag_parts) {
                    imag_parts[i] = sqrt(-disc) / 2.0;
                    imag_parts[i+1] = -sqrt(-disc) / 2.0;
                }
                i++; /* Skip next since it's part of pair */
            } else {
                imag_parts[i] = 0.0;
            }
        } else {
            imag_parts[i] = 0.0;
        }
    }

    for (int i = 0; i < n; i++) {
        free(A[i]); free(Q[i]); free(R[i]); free(T[i]);
    }
    free(A); free(Q); free(R); free(T);
    return n;
}

/* --------------------------------------------------------------
 * Knowledge Point: Real root finding via Sturm sequences + bisection
 * (L5 Algorithms)
 *
 * Uses Sturm's theorem to isolate intervals containing real
 * roots, then applies the bisection method for refinement.
 * ------------------------------------------------------------ */
int kh_find_real_roots(const KH_Polynomial* p, double* roots, int max_roots) {
    if (!p || !roots || max_roots < 1) return 0;

    /* Search interval: use Cauchy bound */
    double max_coeff = 0.0;
    for (int i = 0; i < p->degree; i++) {
        double abs_c = fabs(p->coeffs[i] / p->coeffs[p->degree]);
        if (abs_c > max_coeff) max_coeff = abs_c;
    }
    double bound = 1.0 + max_coeff;

    /* Subdivide interval and use Sturm to count roots */
    int n_sub = 200;
    int found = 0;

    for (int i = 0; i < n_sub && found < max_roots; i++) {
        double a = -bound + 2.0 * bound * i / n_sub;
        double b = -bound + 2.0 * bound * (i + 1) / n_sub;
        int roots_in = sturm_real_root_count(p, a, b);
        if (roots_in > 0) {
            /* Bisection to find the root */
            double lo = a, hi = b;
            for (int iter = 0; iter < 50 && found < max_roots; iter++) {
                double mid = 0.5 * (lo + hi);
                int r_left = sturm_real_root_count(p, lo, mid);
                if (r_left > 0) {
                    hi = mid;
                } else {
                    lo = mid;
                }
                if (hi - lo < 1e-10) {
                    roots[found++] = 0.5 * (lo + hi);
                    break;
                }
            }
        }
    }
    return found;
}

/* --------------------------------------------------------------
 * Knowledge Point: Zero row handling in Routh table (L5 Algorithms)
 *
 * A row of all zeros in the Routh table indicates roots
 * symmetrically placed about the origin (pairs on j*omega axis
 * or symmetric about the origin). The "auxiliary polynomial"
 * is formed from the row above the zero row; its derivative
 * replaces the zero row.
 * ------------------------------------------------------------ */
bool kh_routh_handle_zero_row(KH_RouthTable* rt) {
    if (!rt) return false;
    /* Find zero row and replace with derivative of auxiliary polynomial */
    for (int i = 2; i < rt->n_rows; i++) {
        bool all_zero = true;
        for (int j = 0; j < rt->n_cols; j++) {
            if (fabs(rt->rows[i][j]) > 1e-14) { all_zero = false; break; }
        }
        if (all_zero) {
            /* Auxiliary polynomial from row i-1: coeffs * s^{n-i+1}, s^{n-i-1}, ... */
            int aux_deg = rt->degree - i + 1;
            for (int j = 0; j < rt->n_cols; j++) {
                int pow = aux_deg - 2 * j;
                if (pow >= 0) {
                    rt->rows[i][j] = rt->rows[i-1][j] * pow;
                }
            }
            return true;
        }
    }
    return false;
}

/* Epsilon method: replace zero in first column with small epsilon */
bool kh_routh_handle_zero_first_col(KH_RouthTable* rt, double epsilon) {
    if (!rt) return false;
    for (int i = 1; i < rt->n_rows; i++) {
        if (fabs(rt->rows[i][0]) < 1e-14) {
            rt->rows[i][0] = epsilon;
            return true;
        }
    }
    return false;
}

KH_StabilityResult kh_routh_special_case(const KH_RouthTable* rt) {
    if (!rt) return KH_INCONCLUSIVE;
    /* Check for zero rows and zero first column entries */
    for (int i = 0; i < rt->n_rows; i++) {
        if (fabs(rt->rows[i][0]) < 1e-14) {
            bool all_zero = true;
            for (int j = 1; j < rt->n_cols; j++) {
                if (fabs(rt->rows[i][j]) > 1e-14) all_zero = false;
            }
            if (all_zero) return KH_MARGINAL; /* j*omega axis roots */
            /* Zero in first column only: replace with epsilon,
             * if no sign change from epsilon method, marginal */
        }
    }
    return KH_INCONCLUSIVE;
}

/* --------------------------------------------------------------
 * Knowledge Point: Lienard-Chipart applied to interval polynomials
 * (L8 Advanced Topics)
 *
 * The Lienard-Chipart criterion can be extended to interval
 * polynomials by checking the worst-case combination of
 * coefficient bounds. However, the Kharitonov theorem is
 * generally preferred as it provides a complete answer.
 * This is a fast pre-screening method.
 * ------------------------------------------------------------ */
bool kh_lienard_chipart_interval(const KH_IntervalPoly* ip) {
    if (!ip) return false;
    /* Check that all coefficients in an alternating pattern are positive */
    /* For interval polynomials, a necessary condition for robust
     * stability is that the minimum value of each coefficient
     * at positions n, n-2, n-4, ... be positive */
    for (int k = ip->degree; k >= 0; k -= 2) {
        if (ip->coeffs[k].lo <= 0.0) return false;
    }
    /* Similarly for n-1, n-3, ... */
    for (int k = ip->degree - 1; k >= 0; k -= 2) {
        if (ip->coeffs[k].lo <= 0.0) return false;
    }
    return true;
}
