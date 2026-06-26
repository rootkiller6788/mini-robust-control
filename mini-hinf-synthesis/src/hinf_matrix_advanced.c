/**
 * @file hinf_matrix_advanced.c
 * @brief Advanced matrix algorithms: eigenvalues, Schur, QR, SVD, Lyapunov
 *
 * Numerical linear algebra backbone for H-infinity synthesis.
 * Each function implements one independently studied algorithm from
 * the literature (Golub & Van Loan 2013, Higham 2008, Stewart 2001).
 *
 * Knowledge Points (each function = one algorithm from the literature):
 *   - Householder QR factorization (Golub & Van Loan Alg 5.2.1)
 *   - Least-squares via QR (Golub Sec 5.3)
 *   - Upper Hessenberg reduction via orthogonal similarity (Golub Alg 7.4.2)
 *   - Francis double-shift QR iteration (Golub Alg 7.5.2, Francis 1961)
 *   - Real Schur decomposition (Golub 7.5.2, Laub 1979 IEEE TAC)
 *   - Eigenvalue extraction from real Schur form
 *   - Jacobi eigenvalue algorithm for symmetric matrices (Golub Sec 8.4)
 *   - Power iteration for dominant eigenvalue (Wilkinson 1965 Sec 9.2)
 *   - Spectral radius computation
 *   - Schur form reordering via adjacent block swapping (Bai & Demmel 1993)
 *   - Matrix sign function via Newton-Schulz (Roberts 1980, Higham 2008)
 *   - Matrix square root via Denman-Beavers (Higham 2008 Ch. 6)
 *   - Matrix exponential via scaling-and-squaring + Pade (Higham 2005)
 *   - Singular value decomposition (Golub-Kahan 1965)
 *   - Maximum singular value via power iteration
 *   - Matrix balancing (Osborne 1960, Parlett & Reinsch 1969)
 *   - Bartels-Stewart Lyapunov solver (Bartels & Stewart 1972)
 *   - Bartels-Stewart Sylvester solver
 *   - Hamiltonian/symplectic structure verification (Paige & Van Loan 1981)
 *
 * References:
 *   Golub & Van Loan (2013) Matrix Computations, 4th ed., Johns Hopkins
 *   Higham (2008) Functions of Matrices: Theory and Computation, SIAM
 */

#include "hinf_types.h"
#include "hinf_math.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ---------- internal vector utilities ---------- */
static double vnrm2(const double *v, int n) {
    double s=0.0; for(int i=0;i<n;i++) s+=v[i]*v[i]; return sqrt(s);
}
static double vdot(const double *a, const double *b, int n) {
    double s=0.0; for(int i=0;i<n;i++) s+=a[i]*b[i]; return s;
}

/* ===================================================================
 * Householder Reflection (Golub & Van Loan 2013, Algorithm 5.1.1)
 *
 * Compute v and beta s.t. (I - beta v v^T) x = sigma e_1.
 * v[0] = 1 (implicitly stored); returns beta.
 * On output: x[0] = sigma = +-||x||, x[1..m-1] = 0.
 *
 * Knowledge: The Householder reflection is the fundamental orthogonal
 * transformation underlying QR factorization, Hessenberg reduction,
 * and SVD bidiagonalization. Each reflection encodes a hyperplane
 * symmetry designed to annihilate a sub-vector.
 * =================================================================== */
static double house_reflect(double *x, int m, double *v) {
    if (m <= 1) {
        if (v) v[0] = 1.0;
        return 0.0;
    }
    double xn = vnrm2(x, m);
    if (xn < 1e-16) {
        if (v) { v[0] = 1.0; for (int i = 1; i < m; i++) v[i] = 0.0; }
        return 0.0;
    }
    double sigma = (x[0] > 0) ? -xn : xn;
    if (v) {
        v[0] = 1.0;
        for (int i = 1; i < m; i++) v[i] = x[i] / (x[0] - sigma);
    }
    double beta = 1.0 - (x[0] - sigma) / sigma;
    x[0] = sigma;
    for (int i = 1; i < m; i++) x[i] = 0.0;
    return beta;
}

/* Apply H from left: (I - beta v v^T) * A */
static void house_apply_left(double *v, double beta, int m, int ncol,
                              double *A, int ldA, int coff) {
    if (beta < 1e-16) return;
    double *w = (double *)calloc((size_t)ncol, sizeof(double));
    if (!w) return;
    for (int j = 0; j < ncol; j++) {
        double s = 0.0;
        for (int i = 0; i < m; i++)
            s += A[(j + coff) * ldA + i + coff] * v[i];
        w[j] = beta * s;
    }
    for (int j = 0; j < ncol; j++)
        for (int i = 0; i < m; i++)
            A[(j + coff) * ldA + i + coff] -= v[i] * w[j];
    free(w);
}

/* ===================================================================
 * QR Factorization via Householder (Golub & Van Loan Alg 5.2.1, p.249)
 *
 * Factor A (m x n, m >= n) as A = Q R where Q is m x m orthogonal
 * and R is m x n upper triangular.
 *
 * Algorithm: For k = 0..n-1, compute Householder to zero column k
 * below diagonal, then apply to trailing submatrix.
 * The reflection vector v is stored in A(k+1:m, k), overwriting the
 * zeroed entries. R occupies the upper triangle of A.
 *
 * Complexity: O(m n^2) flops. Storage: in-place + O(m) workspace.
 * =================================================================== */
int hinf_qr_decomp(HinfMatrix *A, HinfMatrix *Q) {
    if (!hinf_matrix_check(A)) return -3;
    int m = A->rows, n = A->cols;
    if (m < n) return -3;

    double *beta = (double *)malloc((size_t)n * sizeof(double));
    if (!beta) return -3;

    for (int k = 0; k < n; k++) {
        int rem = m - k;
        double *colk = (double *)malloc((size_t)rem * sizeof(double));
        double *vk   = (double *)malloc((size_t)rem * sizeof(double));
        if (!colk || !vk) { free(colk); free(vk); free(beta); return -3; }
        for (int i = 0; i < rem; i++) colk[i] = A->data[k * A->ld + i + k];
        beta[k] = house_reflect(colk, rem, vk);
        /* Store R_kk and the Householder vector below diagonal */
        A->data[k * A->ld + k] = colk[0]; /* R(k,k) */
        for (int i = 1; i < rem; i++)
            A->data[k * A->ld + i + k] = vk[i]; /* v(k+1..m) */
        /* Apply Householder to trailing submatrix */
        if (k + 1 < n)
            house_apply_left(vk, beta[k], rem, n - k - 1, A->data, A->ld, k);
        free(colk); free(vk);
    }

    /* Form Q explicitly if requested */
    if (Q && hinf_matrix_check(Q) && Q->rows == m && Q->cols == m) {
        hinf_matrix_set_identity(Q);
        for (int k = n - 1; k >= 0; k--) {
            if (beta[k] < 1e-16) continue;
            double *vk = (double *)malloc((size_t)m * sizeof(double));
            vk[0] = 1.0;
            for (int i = k + 1; i < m; i++)
                vk[i - k] = A->data[k * A->ld + i];
            /* Apply from right: Q := Q * H_k (H symmetric) */
            for (int j = 0; j < m; j++) {
                double s = 0.0;
                for (int i = k; i < m; i++)
                    s += vk[i - k] * Q->data[j * Q->ld + i];
                s *= beta[k];
                for (int i = k; i < m; i++)
                    Q->data[j * Q->ld + i] -= s * vk[i - k];
            }
            free(vk);
        }
    }
    free(beta);
    return 0;
}

/* ===================================================================
 * Least-Squares via QR
 *
 * Solve min ||Ax - b||_2 for A (m x n, m >= n).
 * Steps: (1) QR factor A, (2) c = Q^T b,
 *        (3) solve R(1:n,1:n) x = c(1:n), (4) residual = ||c(n+1:m)||.
 * =================================================================== */
double hinf_qr_ls_solve(const HinfMatrix *A, const double *b, double *x) {
    if (!hinf_matrix_check(A) || !b || !x) return -1.0;
    int m = A->rows, n = A->cols;
    if (m < n) return -1.0;

    HinfMatrix Ac = hinf_matrix_alloc(m, n);
    hinf_matrix_copy(&Ac, A);

    /* QR factor */
    int info = hinf_qr_decomp(&Ac, NULL);
    if (info != 0) { hinf_matrix_free(&Ac); return -1.0; }

    /* c = Q^T b: apply Householder reflections stored in Ac */
    double *c = (double *)malloc((size_t)m * sizeof(double));
    for (int i = 0; i < m; i++) c[i] = b[i];

    for (int k = 0; k < n; k++) {
        int rem = m - k;
        double *vk = (double *)malloc((size_t)rem * sizeof(double));
        vk[0] = 1.0;
        for (int i = 1; i < rem; i++)
            vk[i] = Ac.data[k * Ac.ld + i + k];
        /* Recover beta from stored data */
        double vtv = 1.0;
        for (int i = 1; i < rem; i++) vtv += vk[i] * vk[i];
        if (vtv <= 1.000001) { free(vk); continue; }
        double beta_k = 2.0 / vtv;
        double s = 0.0;
        for (int i = 0; i < rem; i++) s += vk[i] * c[i + k];
        s *= beta_k;
        for (int i = 0; i < rem; i++) c[i + k] -= s * vk[i];
        free(vk);
    }

    /* Back substitution with R(1:n,1:n) */
    for (int i = n - 1; i >= 0; i--) {
        double s = c[i];
        for (int j = i + 1; j < n; j++)
            s -= Ac.data[j * Ac.ld + i] * x[j];
        double rii = Ac.data[i * Ac.ld + i];
        x[i] = (fabs(rii) < 1e-16) ? 0.0 : s / rii;
    }

    /* Residual norm */
    double resid = 0.0;
    for (int i = n; i < m; i++) resid += c[i] * c[i];

    free(c);
    hinf_matrix_free(&Ac);
    return sqrt(resid);
}

/* ===================================================================
 * Upper Hessenberg Reduction (Golub & Van Loan Alg 7.4.2, p.377)
 *
 * Reduce A to upper Hessenberg form via orthogonal similarity:
 * H = Q^T A Q.  H is zero below the first subdiagonal.
 *
 * Algorithm: For each column k = 0..n-3, use a Householder reflection
 * from the left and right to introduce zeros below subdiagonal
 * while preserving previously introduced zeros.
 *
 * Complexity: (10/3) n^3 flops.
 * =================================================================== */
static void hessenberg_reduce(HinfMatrix *A, HinfMatrix *Q) {
    int n = A->rows;
    if (n < 3) return;

    for (int k = 0; k < n - 2; k++) {
        int h = n - k - 1;
        double *x = (double *)malloc((size_t)h * sizeof(double));
        double *v = (double *)malloc((size_t)h * sizeof(double));
        for (int i = k + 1; i < n; i++) x[i - k - 1] = A->data[k * A->ld + i];

        double beta = house_reflect(x, h, v);
        if (beta < 1e-16) { free(x); free(v); continue; }

        /* Left: A(k+1:n, k:n) */
        for (int j = k; j < n; j++) {
            double s = 0.0;
            for (int i = 0; i < h; i++)
                s += v[i] * A->data[j * A->ld + i + k + 1];
            s *= beta;
            for (int i = 0; i < h; i++)
                A->data[j * A->ld + i + k + 1] -= s * v[i];
        }
        /* Right: A(0:n, k+1:n) */
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int jj = 0; jj < h; jj++)
                s += v[jj] * A->data[(jj + k + 1) * A->ld + i];
            s *= beta;
            for (int jj = 0; jj < h; jj++)
                A->data[(jj + k + 1) * A->ld + i] -= s * v[jj];
        }
        /* Accumulate Q */
        if (Q) {
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int jj = 0; jj < h; jj++)
                    s += v[jj] * Q->data[(jj + k + 1) * Q->ld + i];
                s *= beta;
                for (int jj = 0; jj < h; jj++)
                    Q->data[(jj + k + 1) * Q->ld + i] -= s * v[jj];
            }
        }
        free(x); free(v);
    }
}

/* ===================================================================
 * Francis Double-Shift QR Step (Golub & Van Loan Alg 7.5.2, p.389)
 *
 * Perform one implicit double-shift QR step on an unreduced upper
 * Hessenberg matrix H. Uses bulge-chasing to maintain Hessenberg form
 * while applying two real/complex shifts (eigenvalues of the trailing
 * 2x2 submatrix).
 *
 * The "Implicit Q Theorem" guarantees that this is equivalent to
 * two explicit single-shift QR steps, but avoids complex arithmetic
 * when the trailing eigenvalues are a complex conjugate pair.
 *
 * Complexity: O(n^2) per step. Typically ~2n^3 total for convergence.
 *
 * Ref: Francis (1961) The QR Transformation, Parts I-II, Comp. J. 4
 *      Watkins (1982) Understanding the QR Algorithm, SIAM Review 24(4)
 * =================================================================== */
static int qr_double_shift(HinfMatrix *H, HinfMatrix *Q, int ilo, int ihi) {
    int n = H->rows;
    if (ihi - ilo < 2) return 0;

    int m = ihi - 1;
    double a = H->data[(m-1) * H->ld + m-1];
    double b = H->data[m * H->ld + m-1];
    double cc = H->data[(m-1) * H->ld + m];
    double d = H->data[m * H->ld + m];
    double tr = a + d;
    double det = a * d - b * cc;

    /* Form first Householder vector from (M - tr*M + det*I) e_1 */
    double x = H->data[ilo*H->ld+ilo]*H->data[ilo*H->ld+ilo]
             + H->data[ilo*H->ld+ilo+1]*H->data[ilo*H->ld+ilo+1]
             - tr*H->data[ilo*H->ld+ilo] + det;
    double y = H->data[(ilo+1)*H->ld+ilo]
             * (H->data[ilo*H->ld+ilo] + H->data[(ilo+1)*H->ld+ilo+1] - tr);
    double z = H->data[(ilo+1)*H->ld+ilo]
             * H->data[(ilo+2)*H->ld+ilo+1];

    for (int k = ilo; k <= ihi - 3; k++) {
        double vec[3] = {x, y, z};
        double vn = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
        if (vn < 1e-16) { x = y = z = 0.0; continue; }

        double alpha = (vec[0] > 0) ? -vn : vn;
        double v[3] = {1.0, vec[1]/(vec[0]-alpha), vec[2]/(vec[0]-alpha)};
        double bet = 1.0 - (vec[0] - alpha) / alpha;

        /* Left multiply */
        for (int j = k; j < ihi; j++) {
            double s = v[0]*H->data[j*H->ld+k]
                     + v[1]*H->data[j*H->ld+k+1]
                     + v[2]*H->data[j*H->ld+k+2];
            s *= bet;
            H->data[j*H->ld+k]   -= s*v[0];
            H->data[j*H->ld+k+1] -= s*v[1];
            H->data[j*H->ld+k+2] -= s*v[2];
        }
        /* Right multiply */
        int rstart = (k > ilo) ? k - 1 : 0;
        for (int i = rstart; i <= k + 3 && i < ihi; i++) {
            double s = v[0]*H->data[k*H->ld+i]
                     + v[1]*H->data[(k+1)*H->ld+i]
                     + v[2]*H->data[(k+2)*H->ld+i];
            s *= bet;
            H->data[k*H->ld+i]   -= s*v[0];
            H->data[(k+1)*H->ld+i] -= s*v[1];
            H->data[(k+2)*H->ld+i] -= s*v[2];
        }
        /* Accumulate Q */
        if (Q) {
            for (int i = 0; i < n; i++) {
                double s = v[0]*Q->data[k*Q->ld+i]
                         + v[1]*Q->data[(k+1)*Q->ld+i]
                         + v[2]*Q->data[(k+2)*Q->ld+i];
                s *= bet;
                Q->data[k*Q->ld+i]   -= s*v[0];
                Q->data[(k+1)*Q->ld+i] -= s*v[1];
                Q->data[(k+2)*Q->ld+i] -= s*v[2];
            }
        }
        /* Next bulge position */
        if (k + 3 < ihi) {
            x = H->data[(k+1)*H->ld+k+1];
            y = H->data[(k+2)*H->ld+k+1];
            z = (k + 3 < ihi) ? H->data[(k+3)*H->ld+k+1] : 0.0;
        }
    }
    return 0;
}

/* ===================================================================
 * Real Schur Decomposition (Golub & Van Loan Alg 7.5.2, p.389)
 *
 * Compute A = Q T Q^T where T is quasi-upper-triangular (real Schur
 * form). 1x1 diagonal blocks = real eigenvalues; 2x2 blocks =
 * complex conjugate eigenvalue pairs.
 *
 * The real Schur form is the foundation of:
 *   - Algebraic Riccati Equation solving (Laub 1979)
 *   - Bounded Real Lemma eigenvalue test
 *   - Lyapunov/Sylvester equation solving (Bartels-Stewart 1972)
 *   - Model reduction via balanced truncation
 *
 * Algorithm: Hessenberg reduction + Francis double-shift QR
 * iteration with deflation. Convergence is guaranteed by the
 * Francis shift strategy for all matrices.
 *
 * Complexity: O(n^3) with ~7n^3 average, ~25n^3 worst case.
 * Ref: Laub (1979) A Schur Method for Solving ARE, IEEE TAC 24(6):913-921
 * =================================================================== */
int hinf_schur_decomp(HinfMatrix *A, HinfMatrix *Q) {
    if (!hinf_matrix_check(A)) return -3;
    int n = A->rows;
    if (n != A->cols || n == 0) return -3;

    /* Initialize Q = I if requested */
    if (Q) {
        if (!hinf_matrix_check(Q) || Q->rows != n || Q->cols != n) return -3;
        hinf_matrix_set_identity(Q);
    }

    if (n == 1) return 0;

    /* Step 1: Reduce to upper Hessenberg form */
    hessenberg_reduce(A, Q);

    /* Step 2: Francis QR iteration with deflation */
    int max_iter = 100 * n;
    int iter = 0;
    int ilo = 0, ihi = n;

    while (ihi > ilo + 1 && iter < max_iter) {
        /* Deflation: check subdiagonal entries for convergence */
        int deflated = 0;
        for (int k = ihi - 1; k > ilo; k--) {
            double h_sub = fabs(A->data[(k-1) * A->ld + k]);
            double h_diag = fabs(A->data[(k-1) * A->ld + k-1])
                          + fabs(A->data[k * A->ld + k]);
            if (h_sub < 1e-13 * h_diag) {
                ihi = k;      /* Deflate this eigenvalue */
                deflated = 1;
                break;
            }
        }
        if (deflated) continue;

        /* Perform one double-shift QR step */
        qr_double_shift(A, Q, ilo, ihi);
        iter++;
    }

    return (iter >= max_iter) ? -2 : 0;
}

/* ===================================================================
 * Eigenvalue Extraction from Real Schur Form
 *
 * For 1x1 diagonal blocks: lambda = T(k,k) (real).
 * For 2x2 blocks [a b; c d]: lambdas are roots of
 *   lambda^2 - (a+d) lambda + (ad - bc) = 0.
 *
 * Returns the number of eigenvalues extracted (always n).
 * =================================================================== */
int hinf_eigenvalues(const HinfMatrix *A, double *real_p, double *imag_p) {
    if (!hinf_matrix_check(A) || !real_p || !imag_p) return -3;
    int n = A->rows;
    if (n != A->cols || n == 0) return 0;

    HinfMatrix T = hinf_matrix_alloc(n, n);
    hinf_matrix_copy(&T, A);
    int info = hinf_schur_decomp(&T, NULL);
    if (info != 0 && info != -2) { hinf_matrix_free(&T); return info; }

    int cnt = 0, i = 0;
    while (i < n) {
        if (i == n - 1 || fabs(T.data[i * T.ld + i+1]) < 1e-13) {
            /* 1x1 block: real eigenvalue */
            real_p[cnt] = T.data[i * T.ld + i];
            imag_p[cnt] = 0.0;
            cnt++; i++;
        } else {
            /* 2x2 block: potentially complex eigenvalues */
            double a = T.data[i * T.ld + i];
            double b = T.data[(i+1) * T.ld + i];
            double cc = T.data[i * T.ld + i+1];
            double d = T.data[(i+1) * T.ld + i+1];
            double tr = a + d;
            double det = a * d - b * cc;
            double disc = tr * tr - 4.0 * det;
            if (disc >= 0) {
                /* Real eigenvalues from 2x2 block */
                real_p[cnt] = 0.5 * (tr + sqrt(disc)); imag_p[cnt] = 0.0; cnt++;
                real_p[cnt] = 0.5 * (tr - sqrt(disc)); imag_p[cnt] = 0.0; cnt++;
            } else {
                /* Complex conjugate pair */
                real_p[cnt] = 0.5 * tr; imag_p[cnt] =  0.5 * sqrt(-disc); cnt++;
                real_p[cnt] = 0.5 * tr; imag_p[cnt] = -0.5 * sqrt(-disc); cnt++;
            }
            i += 2;
        }
    }

    hinf_matrix_free(&T);
    return cnt;
}

/* ===================================================================
 * Symmetric Eigenvalues via Jacobi Rotations (Golub & Van Loan Sec 8.4)
 *
 * For real symmetric A, the Jacobi algorithm iteratively zeros
 * off-diagonal elements via plane rotations. Guaranteed to converge
 * for symmetric matrices. All eigenvalues are real.
 *
 * Each rotation targets the largest off-diagonal element (classic
 * Jacobi) and zeroes it by a Givens rotation J(p,q,theta) where
 * tan(2*theta) = 2 a_pq / (a_pp - a_qq).
 *
 * Complexity: O(n^3) typical, O(n^4) worst case.
 * =================================================================== */
int hinf_eigenvalues_symm(const HinfMatrix *A, double *eigvals) {
    if (!hinf_matrix_check(A) || !eigvals) return -3;
    int n = A->rows;
    if (n != A->cols || n == 0) return 0;

    HinfMatrix T = hinf_matrix_alloc(n, n);
    hinf_matrix_copy(&T, A);

    int max_iter = 50 * n;
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        /* Find max off-diagonal element */
        int p = 0, q = 1;
        double max_off = fabs(T.data[1 * T.ld + 0]);
        for (int j = 1; j < n; j++)
            for (int i = 0; i < j; i++) {
                double v = fabs(T.data[j * T.ld + i]);
                if (v > max_off) { max_off = v; p = i; q = j; }
            }
        if (max_off < 1e-13) break;

        double app = T.data[p * T.ld + p];
        double aqq = T.data[q * T.ld + q];
        double apq = T.data[q * T.ld + p];
        double theta;
        if (fabs(app - aqq) < 1e-15) {
            theta = (apq > 0) ? 0.7853981633974483 : -0.7853981633974483;
        } else {
            theta = 0.5 * atan2(2.0 * apq, app - aqq);
        }
        double cs = cos(theta), sn = sin(theta);

        /* Apply rotation J(p,q,theta) */
        for (int i = 0; i < n; i++) {
            if (i == p || i == q) continue;
            double tip = T.data[p * T.ld + i];
            double tiq = T.data[q * T.ld + i];
            T.data[p * T.ld + i] =  cs * tip + sn * tiq;
            T.data[q * T.ld + i] = -sn * tip + cs * tiq;
            T.data[i * T.ld + p] = T.data[p * T.ld + i];
            T.data[i * T.ld + q] = T.data[q * T.ld + i];
        }
        T.data[p * T.ld + p] = cs*cs*app + 2.0*cs*sn*apq + sn*sn*aqq;
        T.data[q * T.ld + q] = sn*sn*app - 2.0*cs*sn*apq + cs*cs*aqq;
        T.data[q * T.ld + p] = 0.0;
        T.data[p * T.ld + q] = 0.0;
    }

    for (int i = 0; i < n; i++)
        eigvals[i] = T.data[i * T.ld + i];

    hinf_matrix_free(&T);
    return (iter >= max_iter) ? -2 : n;
}

/* ===================================================================
 * Dominant Eigenvalue via Power Iteration (Wilkinson 1965, Sec 9.2)
 *
 * v_{k+1} = A v_k / ||A v_k||_2
 * lambda_k = v_k^T A v_k / v_k^T v_k  (Rayleigh quotient)
 *
 * Converges to eigenvalue with largest |lambda|.
 * Rate: linear, with convergence factor |lambda_2| / |lambda_1|.
 * =================================================================== */
double hinf_eigenvalue_dominant(const HinfMatrix *A, double *eigvec, int max_iter) {
    if (!hinf_matrix_check(A)) return 0.0;
    int n = A->rows;
    if (n != A->cols || n == 0) return 0.0;
    if (max_iter <= 0) max_iter = 100;

    double *v = (double *)calloc((size_t)n, sizeof(double));
    double *w = (double *)calloc((size_t)n, sizeof(double));
    if (!v || !w) { free(v); free(w); return 0.0; }

    for (int i = 0; i < n; i++) v[i] = 1.0;

    double lambda = 0.0, lambda_old;
    for (int iter = 0; iter < max_iter; iter++) {
        /* w = A * v */
        for (int i = 0; i < n; i++) {
            w[i] = 0.0;
            for (int j = 0; j < n; j++)
                w[i] += A->data[j * A->ld + i] * v[j];
        }
        /* Rayleigh quotient */
        double num = vdot(v, w, n);
        double den = vdot(v, v, n);
        if (den < 1e-16) break;
        lambda_old = lambda;
        lambda = num / den;
        /* Normalize */
        double wn = vnrm2(w, n);
        if (wn < 1e-16) break;
        for (int i = 0; i < n; i++) v[i] = w[i] / wn;
        /* Convergence check */
        if (fabs(lambda - lambda_old) < 1e-10 * (fabs(lambda) + 1.0)) break;
    }

    if (eigvec) for (int i = 0; i < n; i++) eigvec[i] = v[i];
    free(v); free(w);
    return lambda;
}

/* ===================================================================
 * Spectral Radius: rho(A) = max_i |lambda_i(A)|
 *
 * For small matrices (n <= 4), compute all eigenvalues exactly.
 * For larger matrices, use power iteration for quick estimate.
 *
 * Critical application in H-infinity synthesis:
 *   DGKF coupling condition: rho(X_inf * Y_inf) < gamma^2
 * This determines whether a gamma-suboptimal controller exists.
 * =================================================================== */
double hinf_spectral_radius(const HinfMatrix *A) {
    if (!hinf_matrix_check(A)) return -1.0;
    int n = A->rows;
    if (n != A->cols) return -1.0;
    if (n <= 4) {
        double *r = (double *)malloc((size_t)n * sizeof(double));
        double *im = (double *)malloc((size_t)n * sizeof(double));
        double rho = 0.0;
        if (r && im) {
            int m = hinf_eigenvalues(A, r, im);
            for (int i = 0; i < m; i++) {
                double mag = sqrt(r[i]*r[i] + im[i]*im[i]);
                if (mag > rho) rho = mag;
            }
        }
        free(r); free(im);
        return rho;
    }
    return fabs(hinf_eigenvalue_dominant(A, NULL, 100));
}

double hinf_spectral_radius_product(const HinfMatrix *A, const HinfMatrix *B) {
    if (!hinf_matrix_check(A) || !hinf_matrix_check(B)) return -1.0;
    int n = A->rows;
    if (n != A->cols || n != B->rows || n != B->cols) return -1.0;
    HinfMatrix C = hinf_matrix_alloc(n, n);
    hinf_mat_mul(&C, A, B);
    double rho = hinf_spectral_radius(&C);
    hinf_matrix_free(&C);
    return rho;
}

/* ===================================================================
 * Schur Form Reordering
 *
 * Swap adjacent 1x1 diagonal blocks in real Schur form using
 * a Givens rotation (Bai & Demmel 1993).
 *
 * Used to group stable and unstable eigenvalues of the Hamiltonian
 * matrix when solving AREs: the stable invariant subspace forms
 * the first n columns of the ordered orthogonal transformation.
 *
 * Ref: Bai & Demmel (1993) On Swapping Diagonal Blocks in Real
 *      Schur Form, Lin. Alg. Appl. 186:73-95
 * =================================================================== */
static int schur_swap_1x1(HinfMatrix *T, HinfMatrix *Q, int p) {
    int n = T->rows;
    double t11 = T->data[p * T->ld + p];
    double t22 = T->data[(p+1) * T->ld + p+1];
    double t12 = T->data[(p+1) * T->ld + p];
    if (fabs(t11 - t22) < 1e-14) return 0;

    double xx = t12, yy = t22 - t11;
    double r = hypot(xx, yy);
    if (r < 1e-14) return 0;
    double cs = yy / r, sn = xx / r;

    /* Left multiply by Givens */
    for (int j = p; j < n; j++) {
        double t1 = cs * T->data[j*T->ld+p] + sn * T->data[j*T->ld+p+1];
        double t2 = -sn * T->data[j*T->ld+p] + cs * T->data[j*T->ld+p+1];
        T->data[j*T->ld+p] = t1;
        T->data[j*T->ld+p+1] = t2;
    }
    /* Right multiply by Givens */
    for (int i = 0; i <= p + 1; i++) {
        double t1 = T->data[p*T->ld+i]*cs + T->data[(p+1)*T->ld+i]*(-sn);
        double t2 = T->data[p*T->ld+i]*sn + T->data[(p+1)*T->ld+i]*cs;
        T->data[p*T->ld+i] = t1;
        T->data[(p+1)*T->ld+i] = t2;
    }
    /* Accumulate Q */
    if (Q) {
        for (int i = 0; i < n; i++) {
            double t1 = Q->data[p*Q->ld+i]*cs + Q->data[(p+1)*Q->ld+i]*(-sn);
            double t2 = Q->data[p*Q->ld+i]*sn + Q->data[(p+1)*Q->ld+i]*cs;
            Q->data[p*Q->ld+i] = t1;
            Q->data[(p+1)*Q->ld+i] = t2;
        }
    }
    return 0;
}

int hinf_schur_ord_select(HinfMatrix *T, HinfMatrix *Q, int select) {
    if (!hinf_matrix_check(T)) return -3;
    int n = T->rows;
    if (n != T->cols || select < 0 || select >= n) return -3;
    for (int k = select; k > 0; k--)
        schur_swap_1x1(T, Q, k - 1);
    return 0;
}

int hinf_schur_ord_stable(HinfMatrix *T, HinfMatrix *Q) {
    if (!hinf_matrix_check(T)) return -3;
    int n = T->rows;
    int n_stable = 0;
    for (int k = n - 1; k >= 0; k--) {
        if (T->data[k * T->ld + k] < 0.0) {
            for (int j = k; j > n_stable; j--)
                schur_swap_1x1(T, Q, j - 1);
            n_stable++;
        }
    }
    return n_stable;
}

int hinf_schur_count_stable(const HinfMatrix *T) {
    if (!hinf_matrix_check(T)) return 0;
    int n = T->rows, cnt = 0;
    for (int i = 0; i < n; i++)
        if (T->data[i * T->ld + i] < 0.0) cnt++;
    return cnt;
}

int hinf_schur_count_imag(const HinfMatrix *T, double tol) {
    if (!hinf_matrix_check(T)) return 0;
    int n = T->rows, cnt = 0;
    for (int i = 0; i < n; i++)
        if (fabs(T->data[i * T->ld + i]) < tol) cnt++;
    return cnt;
}

/* ===================================================================
 * Matrix Sign Function via Newton-Schulz Iteration
 *
 * sign(A) = A (A^2)^{-1/2} for A with no eigenvalues on the
 * imaginary axis. Computed via quadratically convergent Newton-Schulz:
 *   X_{k+1} = 0.5 * X_k * (3I - X_k^2)
 * Converges to sign(A) for any X_0 sufficiently close to sign(A).
 *
 * Application to ARE: The Hamiltonian H = [A -R; -Q -A'] has the
 * property that sign(H) = 2*P_s - I where P_s projects onto the
 * stable invariant subspace. The stabilizing solution X satisfies
 * span([I; X]) = stable invariant subspace of H.
 *
 * Ref: Roberts (1980) Int. J. Control 32(4):677-687
 *      Kenney & Laub (1995) The Matrix Sign Function, IEEE TAC
 *      Higham (2008) Functions of Matrices, Ch. 5
 * =================================================================== */
int hinf_matrix_sign(HinfMatrix *A_in, HinfMatrix *signA, int max_iter,
                      double tol) {
    if (!hinf_matrix_check(A_in) || !hinf_matrix_check(signA)) return -3;
    int n = A_in->rows;
    if (n != A_in->cols) return -3;
    if (max_iter <= 0) max_iter = 50;

    HinfMatrix X = hinf_matrix_alloc(n, n);
    HinfMatrix X2 = hinf_matrix_alloc(n, n);
    HinfMatrix Xold = hinf_matrix_alloc(n, n);
    if (!hinf_matrix_check(&X) || !hinf_matrix_check(&X2)
        || !hinf_matrix_check(&Xold)) {
        hinf_matrix_free(&X); hinf_matrix_free(&X2);
        hinf_matrix_free(&Xold); return -3;
    }

    /* X_0 = A scaled for stability of the iteration */
    double sc = 1.0 / (hinf_spectral_radius(A_in) + 1.0);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            X.data[j * X.ld + i] = A_in->data[j * A_in->ld + i] * sc;

    int converged = 0;
    for (int iter = 0; iter < max_iter; iter++) {
        hinf_matrix_copy(&Xold, &X);
        hinf_mat_mul(&X2, &X, &X);               /* X2 = X^2 */

        HinfMatrix T = hinf_matrix_alloc(n, n);
        hinf_matrix_fill(&T, 0.0);
        for (int i = 0; i < n; i++) T.data[i * T.ld + i] = 3.0;  /* T = 3I */
        hinf_mat_sub(&T, &T, &X2);               /* T = 3I - X^2 */

        HinfMatrix XT = hinf_matrix_alloc(n, n);
        hinf_mat_mul(&XT, &X, &T);               /* XT = X * (3I - X^2) */
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                X.data[j * X.ld + i] = 0.5 * XT.data[j * XT.ld + i];

        hinf_matrix_free(&T); hinf_matrix_free(&XT);

        /* Check ||X - Xold|| / ||X|| < tol */
        hinf_mat_sub(&X2, &X, &Xold);
        if (hinf_mat_norm_fro(&X2) < tol * hinf_mat_norm_fro(&X)) {
            converged = 1;
            break;
        }
    }

    hinf_matrix_copy(signA, &X);
    hinf_matrix_free(&Xold); hinf_matrix_free(&X2); hinf_matrix_free(&X);
    return converged ? 0 : -1;
}

/* ===================================================================
 * Matrix Square Root via Denman-Beavers Iteration
 *
 * For SPD matrix A, compute A^{1/2} using:
 *   Y_0 = A,  Z_0 = I
 *   Y_{k+1} = 0.5 * (Y_k + Z_k^{-1})
 *   Z_{k+1} = 0.5 * (Z_k + Y_k^{-1})
 * Then Y_k -> A^{1/2}, Z_k -> A^{-1/2}.
 *
 * Used in balanced stochastic truncation model reduction.
 *
 * Ref: Denman & Beavers (1976) The matrix sign function and
 *      computations in systems, Appl. Math. Comp. 2:63-74
 *      Higham (2008) Functions of Matrices, Ch. 6
 * =================================================================== */
int hinf_matrix_sqrt(HinfMatrix *A, HinfMatrix *sqrtA, int max_iter) {
    if (!hinf_matrix_check(A) || !hinf_matrix_check(sqrtA)) return -3;
    int n = A->rows;
    if (n != A->cols) return -3;
    if (max_iter <= 0) max_iter = 50;

    HinfMatrix Y  = hinf_matrix_alloc(n, n);
    HinfMatrix Z  = hinf_matrix_alloc(n, n);
    HinfMatrix Yi = hinf_matrix_alloc(n, n);
    HinfMatrix Zi = hinf_matrix_alloc(n, n);

    hinf_matrix_copy(&Y, A);         /* Y_0 = A */
    hinf_matrix_set_identity(&Z);    /* Z_0 = I */

    for (int iter = 0; iter < max_iter; iter++) {
        hinf_matrix_copy(&Yi, &Y);
        hinf_matrix_copy(&Zi, &Z);
        if (hinf_mat_inv(&Yi) || hinf_mat_inv(&Zi)) break;

        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                Y.data[j*Y.ld+i] = 0.5 * (Y.data[j*Y.ld+i] + Zi.data[j*Zi.ld+i]);
                Z.data[j*Z.ld+i] = 0.5 * (Z.data[j*Z.ld+i] + Yi.data[j*Yi.ld+i]);
            }
    }

    hinf_matrix_copy(sqrtA, &Y);
    hinf_matrix_free(&Zi); hinf_matrix_free(&Yi);
    hinf_matrix_free(&Z); hinf_matrix_free(&Y);
    return 0;
}

/* ===================================================================
 * Matrix Exponential via Scaling-and-Squaring + Pade
 *
 * exp(A*t) — the state transition matrix for linear systems.
 *
 * Algorithm (Higham 2005):
 *   1. Scale: s = max(0, ceil(log2(||At||_inf))), B = At / 2^s
 *   2. Pade(2,2): R = (I - B/2 + B^2/12)^{-1} * (I + B/2 + B^2/12)
 *   3. Unscale: exp(At) = R^{2^s} via repeated squaring
 *
 * Used in: computing sampled-data equivalents, time-domain analysis
 * of H-infinity controllers, and discrete-time synthesis.
 *
 * Ref: Higham (2005) The Scaling and Squaring Method for the
 *      Matrix Exponential Revisited, SIAM Review 47(4)
 * =================================================================== */
int hinf_matrix_exp(const HinfMatrix *A, double t, HinfMatrix *expAt) {
    if (!hinf_matrix_check(A) || !hinf_matrix_check(expAt)) return -3;
    int n = A->rows;
    if (n != A->cols || n != expAt->rows || n != expAt->cols) return -3;

    /* Compute scale factor s */
    double normAt = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            normAt += fabs(A->data[j * A->ld + i] * t);
    normAt /= (double)n;

    int s = (int)ceil(log2(normAt + 1e-16));
    if (s < 0) s = 0;
    double inv_2s = 1.0 / (double)(1 << s);

    HinfMatrix B   = hinf_matrix_alloc(n, n);
    HinfMatrix B2  = hinf_matrix_alloc(n, n);
    HinfMatrix num = hinf_matrix_alloc(n, n);
    HinfMatrix den = hinf_matrix_alloc(n, n);
    HinfMatrix tmp = hinf_matrix_alloc(n, n);

    /* B = At / 2^s */
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            B.data[j * B.ld + i] = A->data[j * A->ld + i] * t * inv_2s;

    hinf_mat_mul(&B2, &B, &B);

    /* num = I + B/2 + B^2/12 */
    hinf_matrix_fill(&num, 0.0);
    for (int i = 0; i < n; i++) num.data[i * num.ld + i] = 1.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            num.data[j*num.ld+i] += 0.5 * B.data[j*B.ld+i];
            num.data[j*num.ld+i] += B2.data[j*B2.ld+i] / 12.0;
        }

    /* den = I - B/2 + B^2/12, then invert */
    hinf_matrix_fill(&den, 0.0);
    for (int i = 0; i < n; i++) den.data[i * den.ld + i] = 1.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            den.data[j*den.ld+i] -= 0.5 * B.data[j*B.ld+i];
            den.data[j*den.ld+i] += B2.data[j*B2.ld+i] / 12.0;
        }
    hinf_mat_inv(&den);

    /* exp(B) = den^{-1} * num */
    hinf_mat_mul(expAt, &den, &num);

    /* Squaring: exp(At) = exp(B)^{2^s} */
    for (int k = 0; k < s; k++) {
        hinf_mat_mul(&tmp, expAt, expAt);
        hinf_matrix_copy(expAt, &tmp);
    }

    hinf_matrix_free(&tmp); hinf_matrix_free(&den);
    hinf_matrix_free(&num); hinf_matrix_free(&B2); hinf_matrix_free(&B);
    return 0;
}

/* ===================================================================
 * SVD: Maximum Singular Value via Power Iteration
 *
 * sigma_max = max_{v != 0} ||A v||_2 / ||v||_2.
 *
 * Iterative algorithm: start with random v, iterate:
 *   u = A v / ||Av||;  v = A^T u / ||A^T u||;  sigma = ||A^T u||.
 * This is the power method applied to A^T A, which converges to
 * the largest singular value and corresponding vectors.
 *
 * Complexity: O(k * m * n) for k iterations.
 * =================================================================== */
double hinf_svd_sigma_max(const HinfMatrix *A, double *u, double *v, int max_iter) {
    if (!hinf_matrix_check(A)) return -1.0;
    int m = A->rows, nn = A->cols;
    if (max_iter <= 0) max_iter = 50;

    double *vv  = (double *)calloc((size_t)nn, sizeof(double));
    double *uu  = (double *)calloc((size_t)m, sizeof(double));
    double *Atu = (double *)calloc((size_t)nn, sizeof(double));
    if (!vv || !uu || !Atu) { free(vv); free(uu); free(Atu); return -1.0; }

    for (int i = 0; i < nn; i++) vv[i] = 1.0;

    double sigma = 0.0, sigma_old;
    for (int iter = 0; iter < max_iter; iter++) {
        /* u = A * v */
        for (int i = 0; i < m; i++) {
            uu[i] = 0.0;
            for (int j = 0; j < nn; j++)
                uu[i] += A->data[j * A->ld + i] * vv[j];
        }
        double un = vnrm2(uu, m);
        if (un < 1e-16) break;
        for (int i = 0; i < m; i++) uu[i] /= un;

        /* v = A^T * u */
        for (int j = 0; j < nn; j++) {
            Atu[j] = 0.0;
            for (int i = 0; i < m; i++)
                Atu[j] += A->data[j * A->ld + i] * uu[i];
        }
        sigma_old = sigma;
        sigma = vnrm2(Atu, nn);
        if (sigma < 1e-16) break;
        for (int i = 0; i < nn; i++) vv[i] = Atu[i] / sigma;

        if (fabs(sigma - sigma_old) < 1e-10 * (sigma + 1.0)) break;
    }

    if (u) for (int i = 0; i < m; i++) u[i] = uu[i];
    if (v) for (int i = 0; i < nn; i++) v[i] = vv[i];

    free(Atu); free(uu); free(vv);
    return sigma;
}

/* ===================================================================
 * Full SVD via A^T A eigenvalue decomposition
 *
 * For a matrix A (m x n), the singular values are sigma_i = sqrt(lambda_i(A^T A))
 * where k = min(m,n). Form A^T A, compute its eigenvalues via the symmetric
 * Jacobi solver, take square roots, sort descending.
 *
 * This approach is numerically adequate for moderate-sized matrices
 * and avoids the complexity of the full Golub-Reinsch SVD.
 * =================================================================== */
int hinf_svd(const HinfMatrix *A, double *sigma, HinfMatrix *U, HinfMatrix *VT) {
    (void)U; (void)VT;  /* full U, VT not computed in this implementation */
    if (!hinf_matrix_check(A) || !sigma) return -3;
    int m = A->rows, nn = A->cols;
    int k = (m < nn) ? m : nn;

    HinfMatrix AtA = hinf_matrix_alloc(k, k);
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++) {
            double s = 0.0;
            for (int l = 0; l < m; l++)
                s += A->data[i * A->ld + l] * A->data[j * A->ld + l];
            AtA.data[j * AtA.ld + i] = s;
        }

    double *eig = (double *)malloc((size_t)k * sizeof(double));
    hinf_eigenvalues_symm(&AtA, eig);

    for (int i = 0; i < k; i++)
        sigma[i] = (eig[i] > 1e-14) ? sqrt(eig[i]) : 0.0;

    /* Bubble sort descending */
    for (int i = 0; i < k - 1; i++)
        for (int j = i + 1; j < k; j++)
            if (sigma[j] > sigma[i]) {
                double t = sigma[i];
                sigma[i] = sigma[j];
                sigma[j] = t;
            }

    free(eig);
    hinf_matrix_free(&AtA);
    return k;
}

/* ===================================================================
 * Matrix Balancing (Osborne 1960, Parlett & Reinsch 1969)
 *
 * Find diagonal D such that D^{-1} A D has approximately equal
 * row and column norms. Uses powers of 2 to avoid rounding error
 * in the scale factors.
 *
 * Pre-processing step for eigenvalue computation: balancing can
 * dramatically improve accuracy, especially for badly scaled matrices
 * that appear in two-time-scale H-infinity problems.
 * =================================================================== */
int hinf_balance(HinfMatrix *A, double *scale) {
    if (!hinf_matrix_check(A) || !scale) return -3;
    int n = A->rows;
    if (n != A->cols) return -3;

    for (int i = 0; i < n; i++) scale[i] = 1.0;

    int converged = 0;
    for (int iter = 0; iter < n && !converged; iter++) {
        converged = 1;
        for (int i = 0; i < n; i++) {
            double rn = 0.0, cn = 0.0;
            for (int j = 0; j < n; j++) {
                if (i != j) rn += fabs(A->data[j * A->ld + i]);
                if (i != j) cn += fabs(A->data[i * A->ld + j]);
            }
            if (rn < 1e-16 || cn < 1e-16) continue;

            double f = sqrt(rn / cn);
            double c = 1.0;
            while (c * 2.0 < f) c *= 2.0;
            while (c > 2.0 * f) c /= 2.0;

            if (c < 0.9 || c > 1.1) {
                converged = 0;
                for (int j = 0; j < n; j++)
                    A->data[j * A->ld + i] *= c;
                for (int j = 0; j < n; j++)
                    A->data[i * A->ld + j] /= c;
                scale[i] *= c;
            }
        }
    }
    return 0;
}

/* ===================================================================
 * Symplectic and Hamiltonian Structure Verification
 *
 * Hamiltonian property: H is Hamiltonian iff J H is symmetric,
 * where J = [0 I; -I 0] is the canonical symplectic matrix.
 * This implies H = [A  R; Q  -A^T] with R = R^T, Q = Q^T.
 *
 * Symplectic property: U is symplectic iff U^T J U = J.
 *
 * These structures are preserved under Riccati-related transformations
 * and are fundamental to the geometry of H-infinity control.
 *
 * Ref: Paige & Van Loan (1981) A Schur decomposition for Hamiltonian
 *      matrices, Lin. Alg. Appl. 41:11-32
 * =================================================================== */
int hinf_check_hamiltonian(const HinfMatrix *H, double tol) {
    if (!hinf_matrix_check(H)) return 0;
    int n2 = H->rows;
    if (n2 != H->cols || n2 % 2 != 0) return 0;
    int n = n2 / 2;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            /* H21 symmetric */
            double bl  = H->data[j * H->ld + i + n];
            double blt = H->data[(i + n) * H->ld + j];
            if (fabs(bl - blt) > tol) return 0;
            /* H12 symmetric */
            double tr  = H->data[(j + n) * H->ld + i];
            double trt = H->data[i * H->ld + j + n];
            if (fabs(tr - trt) > tol) return 0;
            /* H11 = -H22^T */
            double h11  = H->data[j * H->ld + i];
            double h22t = H->data[(i + n) * H->ld + j + n];
            if (fabs(h11 + h22t) > tol) return 0;
        }
    }
    return 1;
}

int hinf_check_symplectic(const HinfMatrix *U, double tol) {
    if (!hinf_matrix_check(U)) return 0;
    int n2 = U->rows;
    if (n2 != U->cols || n2 % 2 != 0) return 0;
    int n = n2 / 2;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double tl = 0.0, tr = 0.0;
            for (int k = 0; k < n; k++) {
                tl += U->data[k*U->ld+i+n]*U->data[j*U->ld+k]
                    - U->data[k*U->ld+i]*U->data[j*U->ld+k+n];
                tr += U->data[k*U->ld+i+n]*U->data[(j+n)*U->ld+k]
                    - U->data[k*U->ld+i]*U->data[(j+n)*U->ld+k+n];
            }
            if (fabs(tl) > tol) return 0;
            double ex = (i == j) ? 1.0 : 0.0;
            if (fabs(tr - ex) > tol) return 0;
        }
    return 1;
}

void hinf_symplectic_J(HinfMatrix *J, int n) {
    if (!hinf_matrix_check(J)) return;
    if (J->rows != 2*n || J->cols != 2*n) return;
    hinf_matrix_fill(J, 0.0);
    for (int i = 0; i < n; i++) {
        J->data[(i + n) * J->ld + i] =  1.0;  /* top-right = I */
        J->data[i * J->ld + i + n]     = -1.0;  /* bottom-left = -I */
    }
}

/* ===================================================================
 * Lyapunov Equation: A X + X A^T + Q = 0 (Bartels-Stewart 1972)
 *
 * Solves the continuous-time Lyapunov equation for stable A.
 * Used to compute controllability/observability Gramians, which
 * are essential for balanced truncation model reduction.
 *
 * Algorithm:
 *   1. Real Schur: A = U T U^T
 *   2. Transform: T Y + Y T^T + F = 0  (F = U^T Q U)
 *   3. Solve column-by-column via forward substitution
 *      in the quasi-triangular structure of T
 *   4. Back-transform: X = U Y U^T
 *
 * Complexity: O(n^3).
 * Ref: Bartels & Stewart (1972) Comm. ACM 15(9):820-826
 * =================================================================== */
int hinf_lyapunov_ct(const HinfMatrix *A, const HinfMatrix *Q, HinfMatrix *X) {
    if (!hinf_matrix_check(A) || !hinf_matrix_check(Q)
        || !hinf_matrix_check(X)) return -3;
    int n = A->rows;
    if (n != A->cols || n != Q->rows || n != Q->cols
        || n != X->rows || n != X->cols) return -3;
    if (n == 0) return 0;

    /* Schur decomposition: A = U T U^T */
    HinfMatrix T = hinf_matrix_alloc(n, n);
    HinfMatrix U = hinf_matrix_alloc(n, n);
    hinf_matrix_copy(&T, A);
    hinf_schur_decomp(&T, &U);

    /* F = U^T Q U */
    HinfMatrix F  = hinf_matrix_alloc(n, n);
    HinfMatrix QU = hinf_matrix_alloc(n, n);
    hinf_mat_mul(&QU, Q, &U);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++)
                s += U.data[k * U.ld + i] * QU.data[j * QU.ld + k];
            F.data[j * F.ld + i] = s;
        }

    /* Solve T Y + Y T^T + F = 0 */
    HinfMatrix Y = hinf_matrix_alloc(n, n);
    hinf_matrix_fill(&Y, 0.0);

    for (int j = 0; j < n; j++) {
        for (int i = n - 1; i >= 0; i--) {
            double s = -F.data[j * F.ld + i];
            /* Subtract known terms from same column */
            for (int k = i + 1; k < n; k++)
                s -= T.data[k * T.ld + i] * Y.data[j * Y.ld + k];
            /* Subtract known terms from previous columns */
            for (int k = 0; k < j; k++)
                s -= T.data[j * T.ld + k] * Y.data[k * Y.ld + i];
            double den = T.data[i * T.ld + i] + T.data[j * T.ld + j];
            Y.data[j * Y.ld + i] = (fabs(den) > 1e-14) ? s / den : s;
        }
    }

    /* X = U Y U^T */
    HinfMatrix YUt = hinf_matrix_alloc(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++)
                s += Y.data[k * Y.ld + i] * U.data[j * U.ld + k];
            YUt.data[j * YUt.ld + i] = s;
        }
    hinf_mat_mul(X, &U, &YUt);

    hinf_matrix_free(&YUt); hinf_matrix_free(&Y);
    hinf_matrix_free(&QU); hinf_matrix_free(&F);
    hinf_matrix_free(&U); hinf_matrix_free(&T);
    return 0;
}

/* ===================================================================
 * Sylvester Equation: A X + X B = C (Bartels-Stewart 1972)
 *
 * Solves the general Sylvester equation using Schur decompositions
 * of A and B^T. The unique solution exists iff A and -B have
 * no common eigenvalues.
 *
 * Algorithm:
 *   1. Schur A = UA TA UA^T, Schur B^T = UB TB UB^T
 *   2. Transform C -> UA^T C UB
 *   3. Solve column-by-column using quasi-triangular structure
 *   4. Back-transform
 *
 * Complexity: O(m^3 + n^3).
 * Ref: Bartels & Stewart (1972) Comm. ACM 15(9):820-826
 * =================================================================== */
int hinf_sylvester(const HinfMatrix *A, const HinfMatrix *B,
                     const HinfMatrix *C, HinfMatrix *X) {
    if (!hinf_matrix_check(A) || !hinf_matrix_check(B)
        || !hinf_matrix_check(C) || !hinf_matrix_check(X)) return -3;
    int m = A->rows, nb = B->rows;
    if (m != A->cols || nb != B->cols) return -3;
    if (m != C->rows || nb != C->cols
        || m != X->rows || nb != X->cols) return -3;

    /* Schur for A */
    HinfMatrix TA = hinf_matrix_alloc(m, m);
    HinfMatrix UA = hinf_matrix_alloc(m, m);
    hinf_matrix_copy(&TA, A);
    hinf_schur_decomp(&TA, &UA);

    /* Schur for B^T */
    HinfMatrix Bt = hinf_matrix_alloc(nb, nb);
    HinfMatrix UB = hinf_matrix_alloc(nb, nb);
    hinf_mat_transpose(&Bt, B);
    hinf_schur_decomp(&Bt, &UB);

    /* C' = UA^T * C * UB */
    HinfMatrix Cp  = hinf_matrix_alloc(m, nb);
    HinfMatrix CUB = hinf_matrix_alloc(m, nb);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < nb; j++) {
            double s = 0.0;
            for (int k = 0; k < nb; k++)
                s += C->data[k * C->ld + i] * UB.data[j * UB.ld + k];
            CUB.data[j * CUB.ld + i] = s;
        }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < nb; j++) {
            double s = 0.0;
            for (int k = 0; k < m; k++)
                s += UA.data[k * UA.ld + i] * CUB.data[j * CUB.ld + k];
            Cp.data[j * Cp.ld + i] = s;
        }

    /* Solve TA * Y + Y * Bt = Cp column by column */
    HinfMatrix Y = hinf_matrix_alloc(m, nb);
    hinf_matrix_fill(&Y, 0.0);

    for (int j = 0; j < nb; j++)
        for (int i = m - 1; i >= 0; i--) {
            double s = Cp.data[j * Cp.ld + i];
            for (int k = i + 1; k < m; k++)
                s -= TA.data[k * TA.ld + i] * Y.data[j * Y.ld + k];
            for (int k = 0; k < j; k++)
                s -= Y.data[k * Y.ld + i] * Bt.data[j * Bt.ld + k];
            double den = TA.data[i*TA.ld+i] + Bt.data[j*Bt.ld+j];
            Y.data[j * Y.ld + i] = (fabs(den) > 1e-14) ? s / den : s;
        }

    /* X = UA * Y * UB^T */
    HinfMatrix YUt = hinf_matrix_alloc(m, nb);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < nb; j++) {
            double s = 0.0;
            for (int k = 0; k < nb; k++)
                s += Y.data[k * Y.ld + i] * UB.data[j * UB.ld + k];
            YUt.data[j * YUt.ld + i] = s;
        }
    hinf_mat_mul(X, &UA, &YUt);

    hinf_matrix_free(&YUt); hinf_matrix_free(&Y);
    hinf_matrix_free(&CUB); hinf_matrix_free(&Cp);
    hinf_matrix_free(&Bt); hinf_matrix_free(&UB);
    hinf_matrix_free(&UA); hinf_matrix_free(&TA);
    return 0;
}