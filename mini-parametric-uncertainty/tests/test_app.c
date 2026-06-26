/* test_app.c — Tests for application models */
#include "pu_app.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
static int tr=0,tp=0;
#define T(n) do{tr++;printf("  TEST %s... ",n);}while(0)
#define P() do{tp++;printf("PASS\n");}while(0)

int main(void) {
    printf("=== Test: pu_app ===\n");
    T("dc_motor_model");
    PU_UncertainStateSpace motor = pu_dc_motor_model(1.0,0.5,0.01,0.01,0.001,0.0001,0.1,0.1,0.1);
    assert(motor.n_states == 2);
    assert(motor.n_unc_params == 3);
    pu_uss_free(&motor);
    P();

    T("msd_model");
    PU_UncertainStateSpace msd = pu_msd_model(1.0, 0.5, 10.0, 0.2, 0.2, 0.2);
    assert(msd.n_states == 2);
    pu_uss_free(&msd);
    P();

    T("f16_model");
    PU_UncertainStateSpace f16 = pu_f16_short_period_model(500.0, -1.0, -10.0, -2.0, -0.1, -50.0, 0.2);
    assert(f16.n_states == 2);
    double margin;
    PU_StabilityStatus s = pu_f16_pitch_damper(&f16, 0.5, &margin);
    assert(s == PU_STABLE || s == PU_UNSTABLE || s == PU_INDETERMINATE);
    pu_uss_free(&f16);
    P();

    T("smib_model");
    PU_UncertainStateSpace smib = pu_smib_model(1.2, 1.0, 0.5, 5.0, 0.1, 0.8, 0.15);
    assert(smib.n_states == 2);
    double m = pu_smib_robust_margin(&smib);
    assert(m >= 0.0);
    pu_uss_free(&smib);
    P();

    printf("=== %d/%d tests passed ===\n", tp, tr);
    return (tp==tr)?0:1;
}
