/* test_mu_core.c — Test core u-synthesis data structures and operations */

#include "mu_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

static int tests_run = 0;
#define TEST(name) do { printf("  %s... ", name); fflush(stdout); tests_run++; } while(0)
#define OK() do { printf("OK\n"); } while(0)

int main(void) {
    printf("=== test_mu_core ===\n");

    /* L1: Matrix allocation and element access */
    TEST("matrix create/set/get");
    {
        MUMatrix *m = mu_mat_create(3, 4);
        assert(m != NULL);
        assert(m->rows == 3 && m->cols == 4);
        mu_mat_set(m, 0, 0, mu_complex(1.0, 2.0));
        MUComplex v = mu_mat_get(m, 0, 0);
        assert(fabs(v.re - 1.0) < 1e-12 && fabs(v.im - 2.0) < 1e-12);
        mu_mat_free(m);
        OK();
    }

    /* L1: Real matrix */
    TEST("real matrix create/set/get");
    {
        MURealMatrix *rm = mu_mat_create_real(2, 3);
        assert(rm != NULL);
        mu_real_mat_set(rm, 1, 2, 42.0);
        assert(fabs(mu_real_mat_get(rm, 1, 2) - 42.0) < 1e-12);
        mu_mat_free_real(rm);
        OK();
    }

    /* L2: Matrix identity */
    TEST("matrix identity");
    {
        MUMatrix *I = mu_mat_identity(4);
        assert(I != NULL);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                MUComplex v = mu_mat_get(I, i, j);
                assert(fabs(v.re - (i == j ? 1.0 : 0.0)) < 1e-12);
                assert(fabs(v.im) < 1e-12);
            }
        mu_mat_free(I);
        OK();
    }

    /* L2: Matrix arithmetic */
    TEST("matrix add/sub/mul");
    {
        MUMatrix *A = mu_mat_create(2, 2);
        MUMatrix *B = mu_mat_create(2, 2);
        mu_mat_set(A, 0, 0, mu_complex(1,0)); mu_mat_set(A, 0, 1, mu_complex(2,0));
        mu_mat_set(A, 1, 0, mu_complex(3,0)); mu_mat_set(A, 1, 1, mu_complex(4,0));
        mu_mat_set(B, 0, 0, mu_complex(5,0)); mu_mat_set(B, 0, 1, mu_complex(6,0));
        mu_mat_set(B, 1, 0, mu_complex(7,0)); mu_mat_set(B, 1, 1, mu_complex(8,0));

        MUMatrix *C = mu_mat_add(A, B);
        assert(fabs(mu_mat_get(C, 0, 0).re - 6.0) < 1e-12);
        assert(fabs(mu_mat_get(C, 1, 1).re - 12.0) < 1e-12);
        mu_mat_free(C);

        C = mu_mat_mul(A, B);
        assert(fabs(mu_mat_get(C, 0, 0).re - 19.0) < 1e-12);
        assert(fabs(mu_mat_get(C, 0, 1).re - 22.0) < 1e-12);
        mu_mat_free(C);

        mu_mat_free(A); mu_mat_free(B);
        OK();
    }

    /* L3: Matrix transpose */
    TEST("matrix transpose and conjugate transpose");
    {
        MUMatrix *A = mu_mat_create(2, 3);
        mu_mat_set(A, 0, 0, mu_complex(2,1));
        mu_mat_set(A, 0, 1, mu_complex(3,0));
        MUMatrix *At = mu_mat_transpose(A);
        assert(At->rows == 3 && At->cols == 2);
        assert(fabs(mu_mat_get(At, 0, 0).re - 2.0) < 1e-12);
        mu_mat_free(A); mu_mat_free(At);

        MUMatrix *B = mu_mat_create(2, 2);
        mu_mat_set(B, 0, 0, mu_complex(3,4));
        MUMatrix *Bh = mu_mat_conjugate_transpose(B);
        assert(fabs(mu_mat_get(Bh, 0, 0).im + 4.0) < 1e-12);
        mu_mat_free(B); mu_mat_free(Bh);
        OK();
    }

    /* L3: Matrix inverse */
    TEST("matrix inverse 2x2");
    {
        MUMatrix *A = mu_mat_create(2, 2);
        mu_mat_set(A, 0, 0, mu_complex(4,0)); mu_mat_set(A, 0, 1, mu_complex(7,0));
        mu_mat_set(A, 1, 0, mu_complex(2,0)); mu_mat_set(A, 1, 1, mu_complex(6,0));
        MUMatrix *Ainv = mu_mat_inverse(A);
        assert(Ainv != NULL);
        MUMatrix *Icheck = mu_mat_mul(A, Ainv);
        assert(fabs(mu_mat_get(Icheck, 0, 0).re - 1.0) < 1e-6);
        assert(fabs(mu_mat_get(Icheck, 0, 1).re) < 1e-6);
        assert(fabs(mu_mat_get(Icheck, 1, 0).re) < 1e-6);
        assert(fabs(mu_mat_get(Icheck, 1, 1).re - 1.0) < 1e-6);
        mu_mat_free(A); mu_mat_free(Ainv); mu_mat_free(Icheck);
        OK();
    }

    /* L3: SVD and spectral norm */
    TEST("SVD and spectral norm");
    {
        MUMatrix *A = mu_mat_create(3, 2);
        mu_mat_set(A, 0, 0, mu_complex(3,0)); mu_mat_set(A, 0, 1, mu_complex(0,0));
        mu_mat_set(A, 1, 0, mu_complex(0,0)); mu_mat_set(A, 1, 1, mu_complex(4,0));
        mu_mat_set(A, 2, 0, mu_complex(0,0)); mu_mat_set(A, 2, 1, mu_complex(0,0));
        double sn = mu_mat_spectral_norm(A);
        assert(fabs(sn - 4.0) < 0.01);
        mu_mat_free(A);
        OK();
    }

    /* L3: Frobenius norm */
    TEST("Frobenius norm");
    {
        MUMatrix *A = mu_mat_create(2, 2);
        mu_mat_set(A, 0, 0, mu_complex(3,0)); mu_mat_set(A, 0, 1, mu_complex(4,0));
        mu_mat_set(A, 1, 0, mu_complex(0,0)); mu_mat_set(A, 1, 1, mu_complex(0,0));
        double fn = mu_mat_frobenius_norm(A);
        assert(fabs(fn - 5.0) < 1e-12);
        mu_mat_free(A);
        OK();
    }

    /* L4: Eigenvalue computation */
    TEST("eigenvalue computation 2x2");
    {
        MUMatrix *A = mu_mat_create(2, 2);
        mu_mat_set(A, 0, 0, mu_complex(1,0)); mu_mat_set(A, 0, 1, mu_complex(2,0));
        mu_mat_set(A, 1, 0, mu_complex(2,0)); mu_mat_set(A, 1, 1, mu_complex(1,0));
        MUComplex evals[2];
        int n = mu_mat_eig(A, evals);
        assert(n == 2);
        /* eigenvalues: 3 and -1 */
        double sum_abs = fabs(evals[0].re) + fabs(evals[1].re);
        assert(fabs(sum_abs - 4.0) < 0.5);
        mu_mat_free(A);
        OK();
    }

    /* L4: Stability check */
    TEST("system stability check");
    {
        MUSystem *sys = mu_system_create(2, 1, 1);
        mu_real_mat_set(sys->A, 0, 0, -1.0);
        mu_real_mat_set(sys->A, 1, 1, -2.0);
        mu_real_mat_set(sys->B, 0, 0, 1.0);
        mu_real_mat_set(sys->B, 1, 0, 1.0);
        mu_real_mat_set(sys->C, 0, 0, 1.0);
        mu_real_mat_set(sys->C, 0, 1, 1.0);
        assert(mu_system_is_stable(sys) == true);
        mu_real_mat_set(sys->A, 0, 0, 1.0); /* make unstable */
        assert(mu_system_is_stable(sys) == false);
        mu_system_free(sys);
        OK();
    }

    /* L5: Frequency response */
    TEST("frequency response at DC");
    {
        MUSystem *sys = mu_system_create(1, 1, 1);
        mu_real_mat_set(sys->A, 0, 0, -1.0);
        mu_real_mat_set(sys->B, 0, 0, 2.0);
        mu_real_mat_set(sys->C, 0, 0, 3.0);
        MUMatrix *G = mu_system_freqresp(sys, 0.0);
        assert(G != NULL);
        /* DC gain = C * (-A)^{-1} * B + D = 3 * 1 * 2 + 0 = 6 */
        double gain = mu_mat_get(G, 0, 0).re;
        assert(fabs(gain - 6.0) < 1e-6);
        mu_mat_free(G);
        mu_system_free(sys);
        OK();
    }

    /* L2: Frequency grid creation */
    TEST("frequency grid creation");
    {
        MUFrequencyGrid *g = mu_frequency_grid_create(0.1, 100.0, 50);
        assert(g != NULL);
        assert(g->num_points == 50);
        assert(g->omega[0] >= 0.1 - 1e-12);
        assert(g->omega[49] <= 100.0 + 1e-12);
        mu_frequency_grid_free(g);
        OK();
    }

    /* L1: Uncertainty structure */
    TEST("uncertainty structure creation");
    {
        MUUncertaintyStructure *s = mu_unc_struct_create(3);
        assert(s != NULL);
        assert(s->num_blocks == 3);
        mu_unc_struct_add_block(s, 0, MU_UNC_REAL_SCALAR, 1, 1.0);
        mu_unc_struct_add_block(s, 1, MU_UNC_FULL_BLOCK, 2, 1.0);
        assert(s->blocks[0].type == MU_UNC_REAL_SCALAR);
        assert(s->blocks[1].size == 2);
        mu_unc_struct_free(s);
        OK();
    }

    printf("\nAll %d tests passed.\n", tests_run);
    return 0;
}
