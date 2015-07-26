/*==========================================================
 * RibodynamicsModel.c 
 *
 * Takes an input scalar representing time (t) 
 * a vector representing the state (x)
 * and a vector of model parameters (theta)
 *and outputs the time derivates (dx)
 *
 * The calling syntax is:
 *
 *		dx = RibodynamicsModel(t,x, theta)
 *
 *========================================================*/

#include "mex.h"
#include <math.h>

// function to determine the atc forcing
double atc_input(double dT)
{
    double p1 = 120; //min
    double p2 = 120; //min
    double p = p1+p2;
    double tp = dT-floor(dT/p)*p;
    double u = 0;
    
    if (tp <= p1)
    {
        u = 2.0; //ng/mL
    }  
    else
    {
        u = 0.0; //ng/mL
    }
    return u;  
}

// function to determine the IPTG forcing for now this is just a stub returning 1 always
double iptg_input(double dT)
{
    return 1.0;
}

/* The computational routine */
void Dx(double dT, double* pdX, double* pdTheta, double* pdDx)
{
    /* the list of known parameter values */
    //PLtetO1 and PLlacO1 parameters
    double f_tet = 2535; //dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    double f_lac = 620; // dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    double a_tet = 11;  // nM/min (Lutz and Bujard, NAR, 1997)
    double a_lac = 11; //nM/min (Lutz and Bujard, NAR, 1997)
    //degradation parameters
    double delta_g = 0.0005; //1/min (Andersen et al., Appl. Environ. Microbiol., 1998)
    double matur = 0.132; //1/min (Iizuka et al., Anal. Chem., 2011)
    double vz = 100; //nM/min (rate of enzymatic degradation) (Hersch, PNAS, 2004)
    double Kz = 75; //nM (dissociation constant of enzymatic degradation) (Hersch, PNAS, 2004)
    // other parameters
    double z0 = 9; // experimentally determined
    double copies = 300; //(plasmid copy number)
    double delta_sm = 0; //1/min
    // reading in the parameters we are currently guessing
    double f_srna = pdTheta[0];
    double k_on = pdTheta[1];
    double k_off = pdTheta[2];
    double k_hyb = pdTheta[3];
    double delta_m = pdTheta[4];
    double delta_s = pdTheta[5];
    double delta_c = pdTheta[6];
    double mu = pdTheta[7];
    double beta = pdTheta[8];
    double c = pdTheta[9];

    //determining the forcing
    double dT0 = 0;
    double u = 0;
    double v = 0;
    
    if (dT < dT0)
    {
        u = 0;
        v = 0;
    }
    else
    {
        u = atc_input(dT-dT0); //aTc
        v = iptg_input(dT-dT0); //IPTG
    }
    
    double dTolerance = 0.00001; // a tolerance to make 0 comparisons
    
    double fu;
    if (u < dTolerance)
    {
        fu = f_tet; //aTc
    }    
    else
    {
        fu = 1;
    }
    
    double fv;
    if (v < dTolerance)
    {
        fv = f_lac; //IPTG
    }    
    else
    {
        fv = 1;
    }
    
    // state equations

    pdDx[0] = copies*a_tet/fu - mu*pdX[0] - delta_s*pdX[0] - k_on*pdX[0]*pdX[1] + k_off*pdX[2]; //sRNA

    pdDx[1] = copies*a_lac/fv - mu*pdX[1] - delta_m*pdX[1] - k_on*pdX[0]*pdX[1] + k_off*pdX[2]; //mRNA

    pdDx[2] = k_on*pdX[0]*pdX[1] - k_off*pdX[2] - k_hyb*pdX[2] - mu*pdX[2] - (delta_sm)*pdX[2]; //sRNA:mRNA_intermediate

    pdDx[3] = k_hyb*pdX[2] - mu*pdX[3] - (delta_c)*pdX[3]; //sRNA:mRNA_stable

    pdDx[4] = beta*pdX[1] + f_srna*beta*pdX[3] - matur*pdX[4] - mu*pdX[4] - delta_g*pdX[4] - vz*pdX[4]/(Kz + pdX[4] + c*(pdX[5]-z0)); //GFP non-mature

    pdDx[5] = (1/c)* (matur*pdX[4] - mu*(pdX[5] + delta_g)*c*(pdX[5]-z0) - (vz*c*(pdX[5]-z0))/(Kz + pdX[4] + c*(pdX[5]-z0))); //measured fluoresence
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input time must be a scalar.");
    }
    
    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
    }
     
    /* get the value of the scalar input  */
    double dT = mxGetScalar(prhs[0]);

    /* create a pointer to the real data in the state vector  */
    double* pdX = mxGetPr(prhs[1]);
    
    /* create a pointer to the real data in the parameter vector  */
    double* pdTheta = mxGetPr(prhs[2]);

    /* get dimensions of the state vector */
    mwSize nrows = mxGetN(prhs[1]); 
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nrows,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    double* pdDx = mxGetPr(plhs[0]);

    /* call the computational routine */
    Dx(dT,pdX,pdTheta,pdDx);
}