/* Example: Port-Hamiltonian system construction and passivity verification */
#include "pbc_core.h"
#include "pbc_port_hamiltonian.h"
#include <stdio.h>
int main(void){
    printf("=== Port-Hamiltonian System Example ===\n");
    PBC_PH_System* msd=pbc_ph_mass_spring_damper(1.0,2.0,0.5);
    printf("MSD system: n=%d ports=%d\n",msd->n,msd->m);
    printf("Interconnection J:\n");
    pbc_matrix_print(msd->J,"J");
    printf("Dissipation R:\n");
    pbc_matrix_print(msd->R,"R");
    printf("Input g:\n");
    pbc_matrix_print(msd->g,"g");
    printf("Passive: %s\n",pbc_ph_is_passive(msd)?"YES":"NO");
    double dH[]={1.0,0.5},xdot[2];
    pbc_ph_dynamics(msd,dH,NULL,xdot);
    printf("dx/dt = [%.4f, %.4f]\n",xdot[0],xdot[1]);
    double pb=pbc_ph_power_balance(msd,dH,NULL);
    printf("Power balance: %.4f\n",pb);
    printf("Feedback interconnect example:\n");
    PBC_PH_System* p2=pbc_ph_mass_spring_damper(0.5,1.0,0.1);
    PBC_PH_System* fb=pbc_ph_feedback_interconnect(msd,p2);
    printf("Combined: n=%d, passive=%s\n",fb->n,pbc_ph_is_passive(fb)?"YES":"NO");
    pbc_ph_free(fb);pbc_ph_free(p2);pbc_ph_free(msd);
    return 0;
}
