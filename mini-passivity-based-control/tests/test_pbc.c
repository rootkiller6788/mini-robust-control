#include "pbc_core.h"
#include "pbc_passivity.h"
#include "pbc_energy_shaping.h"
#include "pbc_port_hamiltonian.h"
#include "pbc_ida.h"
#include "pbc_el_systems.h"
#include "pbc_applications.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static int tr=0,tp=0,tf=0;
void run_extra_tests(void);
void run_more_tests(void);

int main(void){
    printf("=== PBC Test Suite ===\n");
    PBC_Matrix* m=pbc_matrix_identity(3);
    assert(m);tr++;tp++;
    pbc_matrix_free(m);
    printf("  matrix ops: PASS\n");
    tr++;
    PBC_Matrix* a=pbc_matrix_zeros(2,2),*b=pbc_matrix_zeros(2,2);
    pbc_matrix_set(a,0,0,1);pbc_matrix_set(a,0,1,2);
    pbc_matrix_set(a,1,0,3);pbc_matrix_set(a,1,1,4);
    pbc_matrix_set(b,0,0,5);pbc_matrix_set(b,0,1,6);
    pbc_matrix_set(b,1,0,7);pbc_matrix_set(b,1,1,8);
    PBC_Matrix* c=pbc_matrix_add(a,b);
    assert(c&&fabs(pbc_matrix_get(c,0,0)-6.0)<1e-12);
    assert(fabs(pbc_matrix_get(c,1,1)-12.0)<1e-12);
    tp++;pbc_matrix_free(a);pbc_matrix_free(b);pbc_matrix_free(c);
    printf("  matrix add: PASS\n");
    tr++;
    PBC_Matrix* aa=pbc_matrix_zeros(2,2),*bb=pbc_matrix_zeros(2,2);
    pbc_matrix_set(aa,0,0,1);pbc_matrix_set(aa,0,1,2);
    pbc_matrix_set(aa,1,0,3);pbc_matrix_set(aa,1,1,4);
    pbc_matrix_set(bb,0,0,5);pbc_matrix_set(bb,0,1,6);
    pbc_matrix_set(bb,1,0,7);pbc_matrix_set(bb,1,1,8);
    PBC_Matrix* d=pbc_matrix_mul(aa,bb);
    assert(d&&fabs(pbc_matrix_get(d,0,0)-19.0)<1e-12);
    assert(fabs(pbc_matrix_get(d,1,1)-50.0)<1e-12);
    tp++;pbc_matrix_free(aa);pbc_matrix_free(bb);pbc_matrix_free(d);
    printf("  matrix mul: PASS\n");
    tr++;
    PBC_Matrix* ss=pbc_matrix_zeros(2,2);
    pbc_matrix_set(ss,0,0,1);pbc_matrix_set(ss,0,1,3);
    pbc_matrix_set(ss,1,0,3);pbc_matrix_set(ss,1,1,2);
    assert(pbc_matrix_is_symmetric(ss,1e-10));
    pbc_matrix_set(ss,1,0,5);assert(!pbc_matrix_is_symmetric(ss,1e-10));
    pbc_matrix_free(ss);
    PBC_Matrix* kk=pbc_matrix_zeros(2,2);
    pbc_matrix_set(kk,0,1,1.0);pbc_matrix_set(kk,1,0,-1.0);
    assert(pbc_matrix_is_skew_symmetric(kk,1e-10));pbc_matrix_free(kk);tp++;
    printf("  sym/skew: PASS\n");
    run_extra_tests();
    run_more_tests();
    printf("\n=== %d/%d passed ===\n",tp,tr);
    return tf>0?1:0;
}

/* Additional tests via appends */
void run_extra_tests(void) {
    printf("--- Storage & LTI ---\n");
    PBC_Matrix* P=pbc_matrix_identity(2);
    PBC_StorageFunction* sf=pbc_storage_quadratic_init(2,P);
    assert(sf);double x[]={3.0,4.0};
    assert(fabs(pbc_storage_eval(sf,x)-12.5)<1e-10);
    pbc_storage_free(sf);pbc_matrix_free(P);tp++;
    printf("  storage: PASS\n");tr++;

    PBC_LTI_System* s1=pbc_lti_first_order(1.0,1.0);
    assert(s1&&s1->n==1);pbc_lti_free(s1);tp++;
    printf("  LTI 1st-order: PASS\n");tr++;

    PBC_LTI_System* s2=pbc_lti_double_integrator();
    assert(s2&&pbc_relative_degree(s2)==2);pbc_lti_free(s2);tp++;
    printf("  LTI double-int: PASS\n");tr++;

    PBC_LTI_System* s3=pbc_lti_mass_spring_damper(1.0,0.5,2.0);
    assert(s3);pbc_lti_free(s3);tp++;
    printf("  LTI MSD: PASS\n");tr++;

    PBC_LTI_System* s4=pbc_lti_dc_motor(1.0,0.01,0.001,0.0001,0.05,0.05);
    assert(s4);pbc_lti_free(s4);tp++;
    printf("  LTI DC motor: PASS\n");tr++;

    printf("--- KYP & Energy Shaping ---\n");
    PBC_LTI_System* sys=pbc_lti_first_order(1.0,1.0);
    PBC_Matrix* PP=pbc_matrix_identity(1);pbc_matrix_set(PP,0,0,1.0);
    PBC_KYP_Matrix* kyp=pbc_kyp_build(sys,PP);
    assert(kyp&&pbc_kyp_is_passive(kyp,1e-10));
    pbc_kyp_free(kyp);pbc_matrix_free(PP);pbc_lti_free(sys);tp++;
    printf("  KYP lemma: PASS\n");tr++;

    PBC_Matrix* Kp=pbc_matrix_identity(1),*Kd=pbc_matrix_identity(1);
    pbc_matrix_set(Kp,0,0,10.0);pbc_matrix_set(Kd,0,0,2.0);
    double qd[]={0.5};PBC_EnergyShapingParams* ep=pbc_es_params_create(1,Kp,Kd,qd);
    assert(ep);double q[]={1.0},qdot[]={0.2},tau[1];
    pbc_es_compute_control(ep,q,qdot,NULL,tau);
    assert(fabs(tau[0]-(-5.4))<1e-10);
    pbc_es_params_free(ep);pbc_matrix_free(Kp);pbc_matrix_free(Kd);tp++;
    printf("  energy shaping: PASS\n");tr++;

    PBC_Matrix* M=pbc_matrix_identity(1),*Kpp=pbc_matrix_identity(1);
    pbc_matrix_set(M,0,0,2.0);pbc_matrix_set(Kpp,0,0,10.0);
    double qa[]={1.5},qda[]={0.5},qda2[]={1.0};
    assert(fabs(pbc_closed_loop_energy(M,Kpp,qa,qda,qda2,1)-1.5)<1e-12);
    assert(pbc_energy_derivative(qda,Kpp,1)<=0.0);
    assert(pbc_lasalle_condition(qda,Kpp,1,1e-10));
    pbc_matrix_free(M);pbc_matrix_free(Kpp);tp++;
    printf("  closed-loop energy: PASS\n");tr++;
}

void run_more_tests(void) {
    printf("--- Gravity & pH ---\n");
    double m1=1,m2=0.5,l1=0.5,l2=0.3,lc1=0.25,lc2=0.15;
    double qg[]={0.0,0.0},g[2];
    pbc_gravity_2link_planar(m1,m2,l1,l2,lc1,lc2,9.81,qg,g);
    assert(fabs(g[0]-5.64075)<1e-5);assert(fabs(g[1]-0.73575)<1e-5);
    tp++;printf("  gravity 2link: PASS\n");tr++;

    PBC_PH_System* msd=pbc_ph_mass_spring_damper(1.0,2.0,0.5);
    assert(msd&&msd->n==2&&pbc_ph_is_passive(msd));
    double dH[]={2.0,0.0},u[]={0.0},xdot[2];
    pbc_ph_dynamics(msd,dH,u,xdot);
    assert(fabs(xdot[0])<1e-12&&fabs(xdot[1]+2.0)<1e-12);
    pbc_ph_free(msd);tp++;
    printf("  pH MSD dynamics: PASS\n");tr++;

    PBC_PH_System* rlc=pbc_ph_rlc_circuit(10.0,0.001,0.0001,1);
    assert(rlc);pbc_ph_free(rlc);tp++;
    printf("  pH RLC: PASS\n");tr++;

    PBC_PH_System* mot=pbc_ph_dc_motor(1.0,0.01,0.001,0.0001,0.05,0.05);
    assert(mot&&mot->m==2);pbc_ph_free(mot);tp++;
    printf("  pH DC motor: PASS\n");tr++;

    PBC_PH_System* p1=pbc_ph_mass_spring_damper(1.0,2.0,0.5);
    PBC_PH_System* p2=pbc_ph_mass_spring_damper(0.5,3.0,0.3);
    PBC_PH_System* fb=pbc_ph_feedback_interconnect(p1,p2);
    assert(fb&&fb->n==4);pbc_ph_free(fb);
    PBC_PH_System* par=pbc_ph_parallel_interconnect(p1,p2);
    assert(par);pbc_ph_free(par);
    pbc_ph_free(p1);pbc_ph_free(p2);tp++;
    printf("  pH interconnect: PASS\n");tr++;

    printf("--- IDA-PBC & Robotics ---\n");
    double u_ida=pbc_ida_1dof_solve(1.0,0.0,1.0,10.0,2.0,0.3,0.1,0.0,0.0);
    assert(fabs(u_ida+3.2)<1e-12);tp++;
    printf("  IDA-PBC 1DOF: PASS\n");tr++;

    PBC_Matrix* Jd=pbc_matrix_identity(2);
    pbc_matrix_set(Jd,0,0,0);pbc_matrix_set(Jd,0,1,1);
    pbc_matrix_set(Jd,1,0,-1);pbc_matrix_set(Jd,1,1,0);
    PBC_Matrix* Rd=pbc_matrix_identity(2);pbc_matrix_set(Rd,0,0,0);
    double xs[]={1.0,0.0};PBC_IDA_Structure* ids=pbc_ida_structure_create(2,Jd,Rd,xs);
    assert(ids);pbc_ida_structure_free(ids);
    pbc_matrix_free(Jd);pbc_matrix_free(Rd);tp++;
    printf("  IDA structure: PASS\n");tr++;

    PBC_RobotModel* rob=pbc_robot_1dof_pendulum(1.0,0.5,0.25,0.1);
    assert(rob);double qr[]={0.0},gr[1];
    pbc_robot_gravity_vector(rob,qr,gr);
    assert(fabs(gr[0]-2.4525)<1e-5);pbc_robot_free(rob);tp++;
    printf("  robot pendulum: PASS\n");tr++;

    PBC_LinkParams ll1={1.0,0.5,0.25,0.05,0.0,1.0};
    PBC_LinkParams ll2={0.5,0.3,0.15,0.02,0.0,1.0};
    PBC_RobotModel* rob2=pbc_robot_2dof_planar(&ll1,&ll2);
    assert(rob2&&rob2->n_links==2);pbc_robot_free(rob2);tp++;
    printf("  robot 2DOF: PASS\n");tr++;

    printf("--- Applications ---\n");
    PBC_DCDC_Converter* conv=pbc_dcdc_create(PBC_CONVERTER_BOOST,0.001,0.0001,10.0,12.0,24.0);
    assert(conv);double duty=pbc_dcdc_control(conv,0.1,0.01);
    assert(duty>=0.0&&duty<=1.0);pbc_dcdc_step(conv,duty,1e-6);pbc_dcdc_free(conv);tp++;
    printf("  DC-DC boost: PASS\n");tr++;

    PBC_DC_Motor* motor=pbc_motor_create(PBC_MOTOR_DC_BRUSHED,1.0,0.01,0.001,0.0001,0.05,0.05,24.0);
    assert(motor);double V=pbc_motor_speed_control(motor,100.0,0.1,0.05);
    assert(fabs(V)<=24.0);pbc_motor_step(motor,V,0.001);
    double V2=pbc_motor_position_control(motor,0.0,2.0,0.1,0.5);
    assert(fabs(V2)<=24.0);pbc_motor_free(motor);tp++;
    printf("  motor PBC: PASS\n");tr++;

    PBC_Quadrotor* quad=pbc_quadrotor_create(1.0,0.01,0.01,0.02,0.25,1e-5,1e-6);
    assert(quad);quad->phi=0.1;quad->theta=-0.05;quad->psi=0.2;
    double Kp[]={2.0,2.0,1.0},Kd[]={0.5,0.5,0.3},tau_q[3];
    pbc_quadrotor_attitude_control(quad,0.0,0.0,0.0,Kp,Kd,tau_q);
    assert(tau_q[0]<0.0&&tau_q[1]>0.0);
    double thrusts[4],Fz=quad->mass*quad->g;
    pbc_quadrotor_mixing(quad,Fz,tau_q,thrusts);
    for(int i=0;i<4;i++)assert(thrusts[i]>=0.0);
    pbc_quadrotor_free(quad);tp++;
    printf("  quadrotor PBC: PASS\n");tr++;

    PBC_Maglev* mag=pbc_maglev_create(1.0,2.0,0.01,0.001,0.01);
    assert(mag);mag->gap=0.012;double Vm=pbc_maglev_control(mag,100.0,20.0);
    pbc_maglev_step(mag,Vm,0.0001);pbc_maglev_free(mag);tp++;
    printf("  maglev PBC: PASS\n");tr++;

    PBC_pH_Process* php=pbc_ph_process_create(10.0,0.5,0.1,0.1,1.8e-5,7.0);
    assert(php);double up=pbc_ph_process_control(php,1.0,0.1);
    assert(up>=0.0&&up<=1.0);pbc_ph_process_step(php,up,0.1);pbc_ph_process_free(php);tp++;
    printf("  pH process: PASS\n");tr++;

    double tmet[]={0.0,0.1,0.2,0.3,0.4,0.5},emet[]={1.0,0.8,0.5,0.2,0.05,0.0};
    assert(pbc_metric_IAE(tmet,emet,6)>0.0);
    assert(pbc_metric_ISE(tmet,emet,6)>0.0);
    assert(pbc_metric_ITAE(tmet,emet,6)>0.0);tp++;
    printf("  metrics: PASS\n");tr++;

    PBC_Matrix* D=pbc_matrix_identity(2);
    pbc_matrix_set(D,0,0,0);pbc_matrix_set(D,0,1,1);
    pbc_matrix_set(D,1,0,-1);pbc_matrix_set(D,1,1,0);
    PBC_DiracStructure* ds=pbc_dirac_create(D,2,2);
    assert(ds&&pbc_dirac_is_valid(ds,1e-10));
    pbc_dirac_free(ds);pbc_matrix_free(D);
    PBC_BG_Element* bg=pbc_bg_element_create(0,PBC_BG_RESISTOR,"Rtest",10.0);
    assert(bg);bg->flow[0]=2.0;pbc_bg_element_constitutive(bg);
    assert(fabs(bg->effort[0]-20.0)<1e-12);pbc_bg_element_free(bg);tp++;
    printf("  Dirac/bond graph: PASS\n");tr++;
}
