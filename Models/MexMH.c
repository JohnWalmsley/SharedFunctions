#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
/* Path to CVode to be changed as required above */


/* Input Arguments
 * - T_IN - Simulation Time and Time Step input in the form [T_0:time_step:T1]
 * - Y_IN - Initial Conditions
 * - P_IN - Parameter values the model is to be simulated with */

#define	T_IN	prhs[0]
#define	Y_IN	prhs[1]
#define P_IN    prhs[2]

/* Output Arguments
 * - YP_OUT is a vector containing the probability of being in the open state
 * at each time point */

#define	M_OUT	plhs[0]
#define	H_OUT	plhs[1]

static int f(realtype t, N_Vector y, N_Vector ydot, void* pr0);

/* Mex Function*/
static void MexMHScaled(N_Vector ydot, realtype t, N_Vector y, double* pr, double* T, double* Y0, int M, double* mout, double* hout) {
    
    int N=(M-1);
    
    /*Work out T0 and Tfinal from T_IN*/
    realtype T0 =T[0];
    
    realtype Tfinal = ((M-1)*((T[1])-(T[0])));
    
    /* Define tolerances to be used*/
    N_Vector abstol = NULL;
    
    abstol = N_VNew_Serial(3);
    
    NV_Ith_S(abstol, 0) = 1e-8;
    NV_Ith_S(abstol, 1) = 1e-8;
    NV_Ith_S(abstol, 2) = 1e-8;
    realtype reltol = 1e-8;
    
    /* Set up CVode*/
    int flag, k;
    N_Vector y0 = NULL;
    void* cvode_mem = NULL;
    y0 = N_VNew_Serial(2);
    NV_Ith_S(y0, 0) = Y0[0];
    NV_Ith_S(y0, 1) = Y0[1];
    realtype NState = 2;
    
    /* Solver options*/
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    
    /* Initialize the memory */
    flag = CVodeSetUserData(cvode_mem, pr);
    
    flag = CVodeInit(cvode_mem, &f, T0, y0);
    
    /* Set tolerances and maximum step*/
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    flag = CVDense(cvode_mem, NState);
    flag = CVodeSetMaxStep(cvode_mem, 0.1);
    
    /* Call CVode test for error and add result to output vector yout */
    for (k = 1; k < N; ++k) {
        double tout = (k*(Tfinal)/N);
        
        
        if (CVode(cvode_mem, tout, y0, &t, CV_NORMAL) < 0) {
            fprintf(stderr, "Error in CVode: %d\n", flag);
        }
        
        /*Probability of being in open state is equal to 1-probability of being in any other state*/
        mout[k] = NV_Ith_S(y0, 0);
        hout[k] = NV_Ith_S(y0, 1);
        
    }
    
    /*Free memory*/
    N_VDestroy_Serial(y0);
    CVodeFree(&cvode_mem);
    
}


static int f(realtype t, N_Vector y, N_Vector ydot, void* pr0) {
    double* pr = pr0;
    
    /* Define voltage protocol to be used*/
    double v;
    /* This shift is needed for simulated protocol to match the protocol recorded in experiment, which is shifted by 0.1ms as compared to the original input protocol. Consequently, each step is held for 0.1ms longer in this version of the protocol as compared to the input.*/
    double shift = 0.1;
    /* The number corresponding to the protocol to be simulated is passed at the front of the parameter values vector*/
    double protocol_number = pr[0];
    
    
    /* made_up_1_shifted*/
    if (protocol_number==1) {
        double C[6] = {54, 26, 10, 0.007/(2*M_PI), 0.037/(2*M_PI), 0.19/(2*M_PI)};
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<6500+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2500-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2500-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2500-shift)));
            
        }
        if (t>=6500+shift && t<7000+shift)
        {
            v=-120;
        }
        if (t>= 7000+shift && t<8000+shift)
        {
            v = -80;
        }
    }
    
    /* made_up_2_shifted*/
    if (protocol_number==2) {
        double C[6] = {45, -32, -14, 0.00357/(2*M_PI), 0.057/(2*M_PI), 0.119/(2*M_PI)};
        
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<6500+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2500-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2500-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2500-shift)));
            
        }
        if (t>=6500+shift && t<7000+shift)
        {
            v=-120;
        }
        if (t>= 7000+shift && t<8000+shift)
        {
            v = -80;
        }
    }
    
    
    /* made_up_3_shifted*/
    if (protocol_number==3) {
        double C[6] = {50, 28, -11, 0.03/(2*M_PI), 0.123/(2*M_PI), 0.199/(2*M_PI)};
        
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<6500+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2500-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2500-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2500-shift)));
            
        }
        if (t>=6500+shift && t<7000+shift)
        {
            v=-120;
        }
        if (t>= 7000+shift && t<8000+shift)
        {
            v = -80;
        }
    }
    
    
    
    /* wang eq prop*/
    if (protocol_number==4) {
        double C[6] = {19.4704157,29.27670088,	-43.51162825,	0.000413039,	0.001402155,	0.010001046};
        
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<6500+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2500-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2500-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2500-shift)));
            
        }
        if (t>=6500+shift && t<7000+shift)
        {
            v=-120;
        }
        if (t>= 7000+shift && t<8000+shift)
        {
            v = -80;
        }
    }
    
    
    /* diveroli eq prop*/
    if (protocol_number==5) {
        double C[6] = {52.2240105,	-18.9443991,	-38.51567768,	0.00093839,	0.001000016,	0.016839351};
        
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<6500+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2500-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2500-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2500-shift)));
            
        }
        if (t>=6500+shift && t<7000+shift)
        {
            v=-120;
        }
        if (t>= 7000+shift && t<8000+shift)
        {
            v = -80;
        }
    }
    
    
    /* max diff*/
    if (protocol_number==6) {
        double C[6] = {53.09022186,	-27.35887966,	10.69044491,	0.000456713,	0.001055862,	0.028866675};
        
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<6500+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2500-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2500-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2500-shift)));
            
        }
        if (t>=6500+shift && t<7000+shift)
        {
            v=-120;
        }
        if (t>= 7000+shift && t<8000+shift)
        {
            v = -80;
        }
    }
    
    
    
    
    /* original_sine*/
    if (protocol_number==15) {
        double C[6] = {57, 18, 28, 0.031, 0.112, 0.08};
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<8000+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2900-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2900-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2900-shift)));
            
        }
        if (t>=8000+shift && t<8500+shift)
        {
            v=-120;
        }
        if (t>= 8500+shift && t<9500+shift)
        {
            v = -80;
        }
    }
    /* equal_proportions*/
    if (protocol_number==16) {
        double C[6] = {51, 71, 17, -8e-5, 0.004, 0.005};
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<8000+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2900-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2900-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2900-shift)));
            
        }
        if (t>=8000+shift && t<8500+shift)
        {
            v=-120;
        }
        if (t>= 8500+shift && t<9500+shift)
        {
            v = -80;
        }
    }
    
    
    /* maz_wang_div_diff*/
    if (protocol_number==17) {
        double C[6] = {25.3, 99.7, 4.2, -0.00007, 0.003, -0.16};
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<8000+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2900-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2900-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2900-shift)));
            
        }
        if (t>=8000+shift && t<8500+shift)
        {
            v=-120;
        }
        if (t>= 8500+shift && t<9500+shift)
        {
            v = -80;
        }
    }
    
    
    
    if (protocol_number==7)
    {
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==8){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==9){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==10){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==11){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==12){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==13){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }     if (protocol_number==14){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    /* Parameters from input parameters - as the first element of the vector corresponds to the protocol number, the parameter values start from the second element*/
    
    double P0 = pr[1];
    double P1 = pr[2];
    double P2 = pr[3];
    double P3 = pr[4];
    double P4 = pr[5];
    double P5 = pr[6];
    double P6 = pr[7];
    double P7 = pr[8];
    
    /*Get variables*/
    
    const double m = NV_Ith_S(y, 0);
    const double h = NV_Ith_S(y, 1);
    
    /* Model equations*/
    
    const double k1 = P0*exp(P1*v);
    const double k2 = P2*exp(-P3*v);
    const double k3 = P4*exp(P5*v);
    const double k4 = P6*exp(-P7*v);
    
    const double m_inf = k1 / ( k1 + k2 );
    const double m_tau =  1 / ( k1 + k2 );
    const double h_inf = k4 / ( k3 + k4 );
    const double h_tau =  1 / ( k3 + k4 );
        
    NV_Ith_S(ydot, 0) = ( m_inf - m) / m_tau;
    NV_Ith_S(ydot, 1) = ( h_inf - h ) / h_tau;
        
    return 0;
}




/* Mex function definition */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
        
{
    N_Vector ydot;
    realtype t;
    N_Vector y;
    int M;
    
    /* Pointer to input variables/parameters*/
    
    double* pr;
    pr = mxGetPr(prhs[2]);
    
    M = mxGetN(T_IN);
    double* T;
    T = mxGetPr(prhs[0]);
    
    double* Y0;
    Y0 = mxGetPr(prhs[1]);
    
    double *mout;
    double *hout;
    int Numvar;
    Numvar = mxGetN(Y_IN);
    
    /* Create a matrix for the return vector of open probabilities */
    M_OUT = mxCreateDoubleMatrix(M-1, 1, mxREAL);
    H_OUT = mxCreateDoubleMatrix(M-1, 1, mxREAL);
    
    /* Assign pointer output */
    mout = mxGetPr(M_OUT);
    hout = mxGetPr(H_OUT);
    
    /* Mex function call */
    MexMHScaled(ydot, t, y, pr, T, Y0, M, mout, hout);
    return;
    
}
