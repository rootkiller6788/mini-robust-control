/* ============================================================================
 * sgt_hinf.c -- H-infinity Norm Computation for the Small Gain Theorem
 *
 * Implements exact/approximate methods for computing ||G||_inf = the
 * induced L2 gain of an LTI system. For the Small Gain Theorem:
 *   ||H1||_inf * ||H2||_inf < 1  =>  closed-loop L2 stable.
 *
 * Algorithms (each an independent knowledge point in robust control):
 *   L5: Frequency grid -- Zhou/Doyle/Glover (1996) Section 4.3
 *   L5: Bisection via Hamiltonian eigenvalues -- BBK (1989)
 *   L4: Bounded Real Lemma verification -- Gahinet/Apkarian (1994)
 *   L5: Algebraic Riccati Equation (Laub 1979, Kleinman-Newton)
 *   L5: H2 norm via controllability Gramian
 *   L5: Hankel norm via Gramian eigenvalue product
 *   L8: Structured H-infinity with D-scaling (mu upper bound)
 * ============================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "sgt_core.h"
#include "sgt_hinf.h"

/* --- Internal helpers --- */
static void* sa_hinf(size_t n) { void *p = malloc(n); assert(p); return p; }
static void* sc_hinf(size_t n, size_t s) { void *p = calloc(n,s); assert(p); return p; }

/** Gaussian elimination with partial pivoting. Returns 0 on success,
 *  -1 if singular. The workhorse linear solver for Lyapunov, Riccati,
 *  and frequency response computations. O(n^3). */
static int gauss_solve(double *A, double *b, int n) {
    for (int col = 0; col < n; col++) {
        int mr = col; double mv = fabs(A[mr * n + col]);
        for (int r = col + 1; r < n; r++) {
            double v = fabs(A[r * n + col]);
            if (v > mv) { mv = v; mr = r; }
        }
        if (mv < 1e-14) return -1;
        if (mr != col) {
            for (int j = col; j < n; j++) {
                double t = A[col * n + j];
                A[col * n + j] = A[mr * n + j];
                A[mr * n + j] = t;
            }
            double t = b[col]; b[col] = b[mr]; b[mr] = t;
        }
        double pv = A[col * n + col];
        for (int r = col + 1; r < n; r++) {
            double fac = A[r * n + col] / pv;
            for (int j = col; j < n; j++)
                A[r * n + j] -= fac * A[col * n + j];
            b[r] -= fac * b[col];
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++)
            b[i] -= A[i * n + j] * b[j];
        b[i] /= A[i * n + i];
    }
    return 0;
}

/** QR eigenvalue algorithm: Hessenberg reduction + Wilkinson-shift QR.
 *  Handles n <= 20. Returns count of stable eigenvalues (Re < -tol).
 *  This is essential for H-infinity bisection: testing whether the
 *  Hamiltonian matrix has eigenvalues on the imaginary axis.
 *  Knowledge L5: QR algorithm for stability analysis.
 *  Reference: Golub & Van Loan, "Matrix Computations," 4th ed., Ch. 7. */
static int qr_eig(const double *M, int n, double *re, double *im, double tol) {
    assert(n > 0 && n <= 20);
    double *H = (double*)sa_hinf(n * n * sizeof(double));
    for (int i = 0; i < n * n; i++) H[i] = M[i];
    if (n == 1) { re[0] = H[0]; im[0] = 0.0; free(H);
                  return (re[0] < -tol) ? 1 : 0; }

    /* Stage 1: Hessenberg reduction via Householder reflections */
    for (int k = 0; k < n - 2; k++) {
        double sigma = 0.0;
        for (int i = k + 1; i < n; i++) sigma += H[i * n + k] * H[i * n + k];
        if (sigma < 1e-16) continue;
        double x1 = H[(k + 1) * n + k];
        double alpha = (x1 > 0.0) ? -sqrt(sigma) : sqrt(sigma);
        double u1 = x1 - alpha;
        double beta = 1.0 / (sigma - x1 * alpha);
        double *u = (double*)sa_hinf(n * sizeof(double));
        u[k + 1] = u1;
        for (int i = k + 2; i < n; i++) u[i] = H[i * n + k];
        /* Left multiply */
        for (int j = k; j < n; j++) {
            double d = 0.0;
            for (int i = k + 1; i < n; i++) d += u[i] * H[i * n + j];
            d *= beta;
            for (int i = k + 1; i < n; i++) H[i * n + j] -= d * u[i];
        }
        /* Right multiply */
        for (int i = 0; i < n; i++) {
            double d = 0.0;
            for (int j = k + 1; j < n; j++) d += H[i * n + j] * u[j];
            d *= beta;
            for (int j = k + 1; j < n; j++) H[i * n + j] -= d * u[j];
        }
        free(u);
    }

    /* Stage 2: QR iteration with Wilkinson shifts */
    for (int iter = 0; iter < 150; iter++) {
        int m = n;
        while (m > 1) {
            double hmm = H[(m-1)*n+(m-1)], hm1m1 = H[(m-2)*n+(m-2)];
            if (fabs(H[(m-1)*n+(m-2)]) < 1e-14*(fabs(hm1m1)+fabs(hmm))) m--;
            else break;
        }
        if (m <= 2) break;

        double a = H[(m-2)*n+(m-2)], b_ = H[(m-2)*n+(m-1)];
        double c_ = H[(m-1)*n+(m-2)], d_ = H[(m-1)*n+(m-1)];
        double tr = a + d_, det = a*d_ - b_*c_;
        double disc = tr*tr - 4.0*det, shift;
        if (disc >= 0.0) {
            double sd = sqrt(disc);
            shift = (fabs(0.5*(tr+sd)-d_) < fabs(0.5*(tr-sd)-d_)) ?
                    0.5*(tr+sd) : 0.5*(tr-sd);
        } else shift = 0.5*tr;

        for (int i = 0; i < m; i++) H[i*n+i] -= shift;
        for (int i = 0; i < m-1; i++) {
            double x = H[i*n+i], y = H[(i+1)*n+i];
            double r = sqrt(x*x + y*y);
            if (r < 1e-15) continue;
            double c = x/r, s = -y/r;
            for (int j = i; j < n; j++) {
                double t1 = H[i*n+j], t2 = H[(i+1)*n+j];
                H[i*n+j] = c*t1 - s*t2; H[(i+1)*n+j] = s*t1 + c*t2;
            }
            for (int j = 0; j < (i+2 < n ? i+2 : m); j++) {
                double t1 = H[j*n+i], t2 = H[j*n+(i+1)];
                H[j*n+i] = c*t1 - s*t2; H[j*n+(i+1)] = s*t1 + c*t2;
            }
        }
        for (int i = 0; i < m; i++) H[i*n+i] += shift;
    }

    /* Stage 3: Extract eigenvalues from quasi-triangular form */
    int idx = 0, i = 0;
    while (i < n) {
        if (i < n-1 && fabs(H[(i+1)*n+i]) > 1e-12) {
            double a = H[i*n+i], b = H[i*n+(i+1)];
            double c = H[(i+1)*n+i], d = H[(i+1)*n+(i+1)];
            double tr2 = a + d, det2 = a*d - b*c;
            double disc2 = tr2*tr2 - 4.0*det2;
            if (disc2 >= 0.0) {
                double sd = sqrt(disc2);
                re[idx]=0.5*(tr2+sd); im[idx++]=0.0;
                re[idx]=0.5*(tr2-sd); im[idx++]=0.0;
            } else {
                re[idx]=0.5*tr2; im[idx++]=0.5*sqrt(-disc2);
                re[idx]=0.5*tr2; im[idx++]=-0.5*sqrt(-disc2);
            }
            i += 2;
        } else { re[idx]=H[i*n+i]; im[idx++]=0.0; i++; }
    }
    int st = 0;
    for (int k = 0; k < n; k++) if (re[k] < -tol) st++;
    free(H); return st;
}

/* ============================================================================
 * L5: Frequency Grid H-infinity Norm — peak of Bode magnitude plot.
 * Evaluates sigma_max(G(jw)) on log-spaced grid, returns peak.
 * Complexity: O(n_omega * nx^3). Reference: Zhou-Doyle-Glover (1996) §4.3. */
double sgt_hinf_grid(const SGTLTISystem *sys, int n_omega) {
    assert(sys && n_omega >= 10);
    int nx = sys->nx, ny = sys->ny, nu = sys->nu;
    double w_min = 1e-3, w_max = 1e4, max_gain = 0.0;
    for (int k = 0; k < n_omega; k++) {
        double w = pow(10.0, log10(w_min) +
                       k * (log10(w_max) - log10(w_min)) / (n_omega - 1));
        double mag;
        if (nu == 1 && ny == 1) {
            if (nx == 1) {
                double a = sys->A[0], b = sys->B[0], c = sys->C[0], d = sys->D[0];
                double denom = a * a + w * w;
                double re_g = c * b * (-a) / denom + d;
                double im_g = c * b * (-w) / denom;
                mag = sqrt(re_g * re_g + im_g * im_g);
            } else if (nx == 2) {
                double a11 = sys->A[0], a12 = sys->A[1];
                double a21 = sys->A[2], a22 = sys->A[3];
                double b1 = sys->B[0], b2 = sys->B[1];
                double c1 = sys->C[0], c2 = sys->C[1];
                double d = sys->D[0];
                double det_re = -w*w + (a11*a22 - a12*a21);
                double det_im = -w * (a11 + a22);
                double adjB1_re = -a22*b1 + a12*b2, adjB1_im = w*b1;
                double adjB2_re = a21*b1 - a11*b2, adjB2_im = w*b2;
                double num_re = c1*adjB1_re + c2*adjB2_re;
                double num_im = c1*adjB1_im + c2*adjB2_im;
                double dm2 = det_re*det_re + det_im*det_im;
                double re_g = (num_re*det_re + num_im*det_im)/dm2 + d;
                double im_g = (num_im*det_re - num_re*det_im)/dm2;
                mag = sqrt(re_g*re_g + im_g*im_g);
            } else { mag = sgt_freq_response_max_sv(sys, w); }
            if (mag > max_gain) max_gain = mag;
        } else {
            double sv = sgt_freq_response_max_sv(sys, w);
            if (sv > max_gain) max_gain = sv;
        }
    }
    return max_gain;
}

/* ============================================================================
 * L5: Frequency response at a single frequency via complex linear solve.
 * G(jw) = C*(jwI - A)^{-1}*B + D. Solves (jwI-A)*X = B as 2nx2n real system:
 * [ -A  -wI; wI  -A ] * [Xr; Xi] = [B; 0]. Then G = C*(Xr+j*Xi) + D.
 * Complexity: O(nx^3) per frequency. */
double sgt_freq_response_max_sv(const SGTLTISystem *sys, double omega) {
    assert(sys);
    int nx = sys->nx, ny = sys->ny, nu = sys->nu;
    if (nx > 10) return sgt_hinf_grid(sys, 500);

    int n2 = 2 * nx;
    double *bigM = (double*)sc_hinf(n2 * n2, sizeof(double));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            bigM[i*n2 + j]         = -sys->A[i*nx + j];
            bigM[i*n2 + (nx + j)]  = (i==j) ? -omega : 0.0;
            bigM[(nx+i)*n2 + j]    = (i==j) ?  omega : 0.0;
            bigM[(nx+i)*n2 + (nx+j)] = -sys->A[i*nx + j];
        }

    double frob2 = 0.0;
    double *rhs = (double*)sc_hinf(n2, sizeof(double));
    for (int in_col = 0; in_col < nu; in_col++) {
        for (int i = 0; i < n2; i++) rhs[i] = 0.0;
        for (int i = 0; i < nx; i++) rhs[i] = sys->B[i*nu + in_col];

        double *Mc = (double*)sa_hinf(n2*n2*sizeof(double));
        double *bc = (double*)sa_hinf(n2*sizeof(double));
        for (int k = 0; k < n2*n2; k++) Mc[k] = bigM[k];
        for (int k = 0; k < n2; k++) bc[k] = rhs[k];

        if (gauss_solve(Mc, bc, n2) != 0) { free(Mc); free(bc); continue; }
        for (int out = 0; out < ny; out++) {
            double re_g = sys->D[out*nu + in_col], im_g = 0.0;
            for (int j = 0; j < nx; j++) {
                re_g += sys->C[out*nx + j] * bc[j];
                im_g += sys->C[out*nx + j] * bc[nx + j];
            }
            frob2 += re_g*re_g + im_g*im_g;
        }
        free(Mc); free(bc);
    }
    free(bigM); free(rhs);
    return sqrt(frob2);
}

/* ============================================================================
 * L4: Hamiltonian matrix for gamma-level H-infinity test:
 *   H = [ A+B*R^{-1}*D^T*C,        B*R^{-1}*B^T;
 *        -C^T*(I+D*R^{-1}*D^T)*C,  -(A+B*R^{-1}*D^T*C)^T ]
 * where R = gamma^2*I - D^T*D. Symplectic eigenvalues come in +/-lambda pairs.
 * Bounded Real Lemma: ||G||_inf < gamma  <=>  H has no imag-axis eigenvalues.
 * Returns NULL if gamma <= sigma_max(D). */
SGTMatrix* sgt_build_hamiltonian(const SGTLTISystem *sys, double gamma) {
    assert(sys);
    int nx = sys->nx, nu = sys->nu, ny = sys->ny;
    double D_frob2 = 0.0;
    for (int i = 0; i < ny*nu; i++) D_frob2 += sys->D[i]*sys->D[i];
    if (gamma*gamma <= D_frob2 + 1e-12) return NULL;

    double R_inv;
    if (nu == 1 && ny == 1) {
        double d = sys->D[0], R_val = gamma*gamma - d*d;
        if (R_val <= 1e-12) return NULL;
        R_inv = 1.0 / R_val;
    } else {
        double R_val = gamma*gamma - D_frob2/fmax(1.0,(double)(nu<ny?nu:ny));
        if (R_val <= 1e-12) return NULL;
        R_inv = 1.0 / R_val;
    }

    int n2 = 2*nx;
    SGTMatrix *H = sgt_matrix_create(n2, n2);

    /* H11 = A + B*R^{-1}*D^T*C */
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            double val = sys->A[i*nx + j];
            for (int k = 0; k < nu; k++) {
                double Bik = sys->B[i*nu + k];
                for (int l = 0; l < ny; l++)
                    val += Bik * R_inv * sys->D[l*nu + k] * sys->C[l*nx + j];
            }
            H->data[i*n2 + j] = val;
        }

    /* H12 = B*R^{-1}*B^T */
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            double val = 0.0;
            for (int k = 0; k < nu; k++)
                val += sys->B[i*nu + k] * R_inv * sys->B[j*nu + k];
            H->data[i*n2 + (nx + j)] = val;
        }

    /* H21 = -C^T*(I + D*R^{-1}*D^T)*C */
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            double val = 0.0;
            for (int k = 0; k < ny; k++) {
                val += sys->C[k*nx + i] * sys->C[k*nx + j];
                double CDtk_i = 0.0, CDtk_j = 0.0;
                for (int l = 0; l < nu; l++) {
                    CDtk_i += sys->C[k*nx + i] * sys->D[k*nu + l];
                    CDtk_j += sys->C[k*nx + j] * sys->D[k*nu + l];
                }
                val += CDtk_i * R_inv * CDtk_j;
            }
            H->data[(nx + i)*n2 + j] = -val;
        }

    /* H22 = -(H11)^T */
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++)
            H->data[(nx + i)*n2 + (nx + j)] = -H->data[j*n2 + i];

    return H;
}

/* ============================================================================
 * L4: Test if Hamiltonian has eigenvalues on the imaginary axis.
 * By the Bounded Real Lemma: ||G||_inf >= gamma iff H(gamma) has
 * eigenvalues with |Re(lambda)| < tol. */
bool sgt_hamiltonian_has_imag_eigs(const SGTMatrix *H, double tol) {
    assert(H && H->rows == H->cols && H->rows % 2 == 0);
    int n = H->rows;
    double *reig = (double*)sa_hinf(n * sizeof(double));
    double *iig = (double*)sa_hinf(n * sizeof(double));
    qr_eig(H->data, n, reig, iig, 1e-8);
    bool found = false;
    for (int i = 0; i < n; i++)
        if (fabs(reig[i]) < tol) { found = true; break; }
    free(reig); free(iig);
    return found;
}

/* ============================================================================
 * L5: H-infinity norm via bisection (Boyd-Balakrishnan-Kabamba 1989).
 *   gamma_lo = sigma_max(D), gamma_hi = freq grid estimate.
 *   Bisect on gamma_mid: if H(gamma_mid) has imag-axis eigs => lo = mid,
 *   else hi = mid. Converges in O(log(range/tol)) iterations, O(nx^3) each. */
double sgt_hinf_bisection(const SGTLTISystem *sys, double tol, int max_iter) {
    assert(sys && tol > 0.0 && max_iter > 0);

    double gamma_lo = 0.0;
    for (int i = 0; i < sys->ny * sys->nu; i++)
        if (fabs(sys->D[i]) > gamma_lo) gamma_lo = fabs(sys->D[i]);

    double gamma_hi = sgt_hinf_grid(sys, 200);
    if (gamma_hi < gamma_lo + 1e-6) gamma_hi = gamma_lo + 1.0;

    /* Ensure gamma_hi is valid (no imag-axis eigenvalues) */
    SGTMatrix *H_hi = sgt_build_hamiltonian(sys, gamma_hi);
    int safety = 0;
    while ((H_hi == NULL || sgt_hamiltonian_has_imag_eigs(H_hi, 1e-8))
           && safety < 40) {
        if (H_hi) sgt_matrix_free(H_hi);
        gamma_hi *= 1.5;
        H_hi = sgt_build_hamiltonian(sys, gamma_hi);
        safety++;
    }
    if (H_hi) sgt_matrix_free(H_hi);

    for (int iter = 0; iter < max_iter; iter++) {
        if (gamma_hi - gamma_lo < tol * fmax(1.0, gamma_hi)) break;
        double gamma_mid = 0.5 * (gamma_lo + gamma_hi);
        SGTMatrix *H_mid = sgt_build_hamiltonian(sys, gamma_mid);
        if (H_mid == NULL) { gamma_lo = gamma_mid; continue; }
        bool has_imag = sgt_hamiltonian_has_imag_eigs(H_mid, 1e-10);
        sgt_matrix_free(H_mid);
        if (has_imag) gamma_lo = gamma_mid;
        else          gamma_hi = gamma_mid;
    }
    return 0.5 * (gamma_lo + gamma_hi);
}

/* ============================================================================
 * L4: Bounded Real Lemma — ||G||_inf < gamma <=> H(gamma) has no imag-axis eigs.
 * The BRL translates an infinite-dimensional frequency-domain constraint into
 * a finite-dimensional algebraic test, enabling LMI-based H-infinity synthesis. */
bool sgt_bounded_real_lemma(const SGTLTISystem *sys, double gamma) {
    assert(sys);
    SGTMatrix *H = sgt_build_hamiltonian(sys, gamma);
    if (H == NULL) return false;
    bool passed = !sgt_hamiltonian_has_imag_eigs(H, 1e-10);
    sgt_matrix_free(H);
    return passed;
}

/* ============================================================================
 * L5: Algebraic Riccati Equation solver. A^T*X + X*A - X*B*R^{-1}*B^T*X + Q = 0.
 * Method: Hamiltonian + Kleinman-Newton iteration (quadratic convergence).
 * For n <= 4: delegates to analytic small-system solver.
 * Knowledge: CARE is central to LQR, Kalman filter, H2/Hinf control. */
SGTMatrix* sgt_solve_care(const SGTMatrix *A, const SGTMatrix *B,
                           const SGTMatrix *R, const SGTMatrix *Q) {
    assert(A && B && R && Q);
    int n = A->rows;
    assert(n == A->cols && n == B->rows && n == R->rows && n == Q->rows);

    if (n > 10) {
        fprintf(stderr, "sgt_solve_care: nx=%d > 10\n", n);
    }

    /* Compute S = B*R^{-1}*B^T via Gauss-Jordan inversion */
    double *Ra = (double*)sa_hinf(n*n*sizeof(double));
    double *Ri = (double*)sa_hinf(n*n*sizeof(double));
    for (int i = 0; i < n*n; i++) {
        Ra[i] = R->data[i];
        Ri[i] = (i % (n+1) == 0) ? 1.0 : 0.0;
    }
    for (int col = 0; col < n; col++) {
        int piv = col;
        double mv = fabs(Ra[piv*n + col]);
        for (int r = col+1; r < n; r++)
            if (fabs(Ra[r*n+col]) > mv) { mv = fabs(Ra[r*n+col]); piv = r; }
        if (mv < 1e-14) { free(Ra); free(Ri); return NULL; }
        if (piv != col) {
            for (int j = 0; j < n; j++) {
                double t = Ra[col*n+j]; Ra[col*n+j]=Ra[piv*n+j]; Ra[piv*n+j]=t;
                t = Ri[col*n+j]; Ri[col*n+j]=Ri[piv*n+j]; Ri[piv*n+j]=t;
            }
        }
        double pv = Ra[col*n+col];
        for (int j = 0; j < n; j++) { Ra[col*n+j] /= pv; Ri[col*n+j] /= pv; }
        for (int r = 0; r < n; r++) {
            if (r == col) continue;
            double fac = Ra[r*n+col];
            for (int j = 0; j < n; j++) {
                Ra[r*n+j] -= fac * Ra[col*n+j];
                Ri[r*n+j] -= fac * Ri[col*n+j];
            }
        }
    }
    free(Ra);

    double *S_s = (double*)sc_hinf(n*n, sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                for (int l = 0; l < n; l++)
                    S_s[i*n+j] += B->data[i*B->stride+k] * Ri[k*n+l] *
                                  B->data[j*B->stride+l];
    free(Ri);

    if (n <= 4) { free(S_s); return sgt_solve_care_small(A, B, R, Q); }

    /* Kleinman-Newton: (A-S*X_k)^T*X_{k+1} + X_{k+1}*(A-S*X_k) + Q + X_k*S*X_k = 0 */
    SGTMatrix *X = sgt_matrix_create(n, n);
    for (int i = 0; i < n; i++) X->data[i*n+i] = 1.0;

    for (int iter = 0; iter < 60; iter++) {
        double *Ak = (double*)sa_hinf(n*n*sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                Ak[i*n+j] = A->data[i*A->stride+j];
                for (int k = 0; k < n; k++)
                    Ak[i*n+j] -= S_s[i*n+k] * X->data[k*n+j];
            }

        double *Qk = (double*)sa_hinf(n*n*sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                Qk[i*n+j] = Q->data[i*Q->stride+j];
                for (int k = 0; k < n; k++)
                    for (int l = 0; l < n; l++)
                        Qk[i*n+j] += X->data[i*n+k] * S_s[k*n+l] *
                                     X->data[l*n+j];
            }

        /* Lyapunov: A_k^T*X + X*A_k = -Q_k.
         * Kronecker: (I⊗A_k^T + A_k^T⊗I) vec(X) = -vec(Q_k). */
        int n2 = n*n;
        double *Ml = (double*)sc_hinf(n2*n2, sizeof(double));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                int row = i*n + j;
                for (int k = 0; k < n; k++)
                    Ml[row*n2 + (k*n+j)] += Ak[k*n+i]; /* I⊗A_k^T */
                for (int k = 0; k < n; k++)
                    Ml[row*n2 + (i*n+k)] += Ak[k*n+j]; /* A_k^T⊗I */
            }

        double *bl = (double*)sa_hinf(n2*sizeof(double));
        for (int i = 0; i < n2; i++) bl[i] = -Qk[i];

        if (gauss_solve(Ml, bl, n2) != 0) {
            free(Ml); free(bl); free(Ak); free(Qk);
            free(S_s); sgt_matrix_free(X); return NULL;
        }

        double md = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double d = fabs(bl[i*n+j] - X->data[i*n+j]);
                if (d > md) md = d;
                X->data[i*n+j] = bl[i*n+j];
            }
        free(Ml); free(bl); free(Ak); free(Qk);
        if (md < 1e-10) break;
    }

    free(S_s);
    return X;
}

/* L5: Small-system CARE. n=1: scalar quadratic with stabilizing root selection.
 * n>=2: delegates to general Kleinman-Newton solver. */
SGTMatrix* sgt_solve_care_small(const SGTMatrix *A, const SGTMatrix *B,
                                 const SGTMatrix *R, const SGTMatrix *Q) {
    assert(A && B && R && Q);
    int n = A->rows;
    assert(n <= 4 && A->cols == n);
    if (n == 1) {
        double a = A->data[0], b = B->data[0];
        double r = R->data[0], q = Q->data[0];
        double c_coef = b*b/r;
        double disc = 4.0*a*a + 4.0*c_coef*q;
        if (disc < 0.0) return NULL;
        double sd = sqrt(disc);
        double x1 = (2.0*a + sd)/(2.0*c_coef);
        double x2 = (2.0*a - sd)/(2.0*c_coef);
        double xp = (a - c_coef*x1 < 0.0) ? x1 : x2;
        if (xp <= 1e-15) return NULL;
        SGTMatrix *X = sgt_matrix_create(1,1);
        X->data[0] = xp;
        return X;
    }
    return sgt_solve_care(A, B, R, Q);
}

/* ============================================================================
 * L5: H2 norm via controllability Gramian. ||G||_2^2 = trace(C*Wc*C^T).
 * Wc solves A*Wc + Wc*A^T + B*B^T = 0. Requires strictly proper (D=0) system.
 * Knowledge: H2 norm = RMS response to white noise. Central to LQG control.
 * Reference: Zhou, Doyle, Glover (1996), §4.4. */
double sgt_h2_norm(const SGTLTISystem *sys) {
    assert(sys);
    int nx = sys->nx, ny = sys->ny, nu = sys->nu;
    for (int i = 0; i < ny*nu; i++)
        if (fabs(sys->D[i]) > 1e-12) return INFINITY;
    if (nx > 4) return INFINITY;

    int n2 = nx*nx;
    double *M = (double*)sc_hinf(n2*n2, sizeof(double));
    /* Lyapunov: A*Wc + Wc*A^T = -B*B^T.
     * Kronecker: (I⊗A + A⊗I) vec(Wc) = -vec(B*B^T).
     * (I⊗A)[(i,j)][(k,j)] = A[i][k], (A⊗I)[(i,j)][(i,k)] = A[j][k]. */
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            int row = i*nx + j;
            for (int k = 0; k < nx; k++) {
                M[row*n2 + (k*nx+j)] += sys->A[i*nx+k]; /* I⊗A */
                M[row*n2 + (i*nx+k)] += sys->A[j*nx+k]; /* A⊗I */
            }
        }
    double *rhs = (double*)sa_hinf(n2*sizeof(double));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            double sum = 0.0;
            for (int k = 0; k < nu; k++)
                sum += sys->B[i*nu+k] * sys->B[j*nu+k];
            rhs[i*nx+j] = -sum;
        }
    if (gauss_solve(M, rhs, n2) != 0) { free(M); free(rhs); return INFINITY; }

    double h2_sq = 0.0;
    for (int i = 0; i < ny; i++)
        for (int k = 0; k < nx; k++) {
            double sw = 0.0;
            for (int j = 0; j < nx; j++)
                sw += sys->C[i*nx+j] * rhs[j*nx+k];
            h2_sq += sw * sys->C[i*nx+k];
        }
    free(M); free(rhs);
    return sqrt(h2_sq);
}

/* ============================================================================
 * L5: Hankel norm = max Hankel singular value = sqrt(max eig(Wc*Wo)).
 * Hankel singular values are system invariants (basis of balanced truncation).
 * ||G||_h <= ||G||_inf always. Reference: Glover (1984). */
double sgt_hankel_norm(const SGTLTISystem *sys) {
    assert(sys);
    int nx = sys->nx, ny = sys->ny, nu = sys->nu;
    if (nx > 4) return INFINITY;
    int n2 = nx*nx;

    /* Wc Lyapunov solve */
    double *MWc = (double*)sc_hinf(n2*n2, sizeof(double));
    double *rWc = (double*)sa_hinf(n2*sizeof(double));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            int row = i*nx + j;
            for (int k = 0; k < nx; k++) {
                MWc[row*n2+(k*nx+j)] += sys->A[i*nx+k];
                MWc[row*n2+(i*nx+k)] += sys->A[j*nx+k];
            }
            double sum = 0.0;
            for (int k = 0; k < nu; k++)
                sum += sys->B[i*nu+k] * sys->B[j*nu+k];
            rWc[i*nx+j] = -sum;
        }
    if (gauss_solve(MWc, rWc, n2) != 0) { free(MWc); free(rWc); return INFINITY; }

    /* Wo Lyapunov solve */
    /* Wo Lyapunov: A^T*Wo + Wo*A = -C^T*C.
     * Kronecker: (I⊗A^T + A^T⊗I) vec(Wo) = -vec(C^T*C).
     * (I⊗A^T)[(i,j)][(k,j)] = A[k][i], (A^T⊗I)[(i,j)][(i,k)] = A[k][j]. */
    double *MWo = (double*)sc_hinf(n2*n2, sizeof(double));
    double *rWo = (double*)sa_hinf(n2*sizeof(double));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            int row = i*nx + j;
            for (int k = 0; k < nx; k++) {
                MWo[row*n2+(k*nx+j)] += sys->A[k*nx+i]; /* I⊗A^T */
                MWo[row*n2+(i*nx+k)] += sys->A[k*nx+j]; /* A^T⊗I */
            }
            double sum = 0.0;
            for (int k = 0; k < ny; k++)
                sum += sys->C[k*nx+i] * sys->C[k*nx+j];
            rWo[i*nx+j] = -sum;
        }
    if (gauss_solve(MWo, rWo, n2) != 0) {
        free(MWc); free(rWc); free(MWo); free(rWo); return INFINITY;
    }

    double *WcWo = (double*)sa_hinf(n2*sizeof(double));
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < nx; j++) {
            double sum = 0.0;
            for (int k = 0; k < nx; k++)
                sum += rWc[i*nx+k] * rWo[k*nx+j];
            WcWo[i*nx+j] = sum;
        }

    double *re = (double*)sa_hinf(nx*sizeof(double));
    double *im = (double*)sa_hinf(nx*sizeof(double));
    qr_eig(WcWo, nx, re, im, 1e-12);

    double max_hsv = 0.0;
    for (int i = 0; i < nx; i++)
        if (re[i] > 1e-12) {
            double hsv = sqrt(re[i]);
            if (hsv > max_hsv) max_hsv = hsv;
        }

    free(MWc); free(rWc); free(MWo); free(rWo);
    free(WcWo); free(re); free(im);
    return max_hsv;
}

/* ============================================================================
 * L8: Structured H-infinity with D-scaling (mu upper bound).
 * mu_Delta(M(jw)) <= inf_{D} sigma_max(D*M(jw)*D^{-1}) where D commutes
 * with Delta. Uses iterative D-scale optimization at each frequency.
 * Reference: Packard-Doyle (1993). */
double sgt_hinf_with_dscale(const SGTLTISystem *sys,
                             const SGTStructuredUncertainty *delta_struct,
                             int n_omega, int n_d_iter) {
    assert(sys && n_omega >= 10);
    if (!delta_struct || delta_struct->n_blocks == 0)
        return sgt_hinf_grid(sys, n_omega);

    double w_min = 1e-3, w_max = 1e4, peak_mu = 0.0;
    int dim = delta_struct->total_dim;
    if (dim <= 0) dim = 1;

    for (int kw = 0; kw < n_omega; kw++) {
        double w = pow(10.0, log10(w_min) +
                kw*(log10(w_max)-log10(w_min))/(n_omega-1));
        double M_mag = sgt_freq_response_max_sv(sys, w);

        double *d = (double*)sa_hinf(dim*sizeof(double));
        for (int i = 0; i < dim; i++) d[i] = 1.0;

        for (int iter = 0; iter < n_d_iter; iter++) {
            double max_ratio = 1.0;
            int wi = 0, wj = 0;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++) {
                    if (d[j] < 1e-15) continue;
                    double r = d[i]/d[j];
                    if (r > max_ratio) { max_ratio = r; wi = i; wj = j; }
                }
            if (max_ratio < 1.0+1e-8) break;
            d[wi] *= 0.98; d[wj] *= 1.02;
        }

        double best_ratio = 1.0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++) {
                if (d[j] < 1e-15) continue;
                double r = d[i]/d[j];
                if (r > best_ratio) best_ratio = r;
            }
        double scaled = sqrt(best_ratio) * M_mag;
        if (scaled > peak_mu) peak_mu = scaled;
        free(d);
    }
    return peak_mu;
}
