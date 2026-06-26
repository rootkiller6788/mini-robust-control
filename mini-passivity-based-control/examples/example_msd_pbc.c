/* Example: Mass-Spring-Damper under PBC energy shaping + damping injection */
#include "pbc_core.h"
#include "pbc_energy_shaping.h"
#include <stdio.h>
#include <stdlib.h>
int main(void){
    printf("=== MSD Energy Shaping PBC ===\n");
    int n=1;
    PBC_Matrix* Kp=pbc_matrix_identity(1);
    PBC_Matrix* Kd=pbc_matrix_identity(1);
    pbc_matrix_set(Kp,0,0,10.0);
    pbc_matrix_set(Kd,0,0,2.0);
    double qd[]={0.5};
    PBC_EnergyShapingParams* p=pbc_es_params_create(1,Kp,Kd,qd);
    double q=1.0,qdot=0.2,tau[1];
    pbc_es_compute_control(p,&q,&qdot,NULL,tau);
    printf("q=%.3f qdot=%.3f qd=%.3f => tau=%.4f\n",q,qdot,qd[0],tau[0]);
    PBC_Matrix* M=pbc_matrix_identity(1);
    pbc_matrix_set(M,0,0,1.0);
    double H_cl=pbc_closed_loop_energy(M,Kp,&q,&qdot,qd,n);
    printf("Closed-loop energy: %.4f\n",H_cl);
    double dH=pbc_energy_derivative(&qdot,Kd,n);
    printf("dH/dt = %.4f (<=0 => stable)\n",dH);
    pbc_es_params_free(p);
    pbc_matrix_free(Kp);pbc_matrix_free(Kd);pbc_matrix_free(M);
    printf("Done.\n");
    return 0;
}
