clc;

%% specify an initial guess for the parameters. 
%The order of parameters in the parameter vector is
%[copies,a_lac,a_tet,f_srna,f_lac,f_tet,k_on,k_off,k_hyb,delta_m,delta_s,delta_sm,delta_c,delta_g,matur,mu,beta,c,z0,vz,Kz]

theta = [
    300,%(plasmid copy number)
    11, % a_lac nM/min (Lutz and Bujard, NAR, 1997)
    11,% a_tetnM/min (Lutz and Bujard, NAR, 1997)
    435, % f_srna dimensionless (fold) TO BE FITTED
    620, % f_lac dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    2535, %dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    1e5, %1/(nM*min) TO BE FITTED
    1e7, %1/min TO BE FITTED 
    0.01, %1/min TO BE FITTED
    0.175, %1/min TO BE FITTED
    0.1, %1/min TO BE FITTED
    0, %1/min
    0.1, %1/min TO BE FITTED
    0.0005, %1/min (Andersen et al., Appl. Environ. Microbiol., 1998)
    0.132, %1/min (Iizuka et al., Anal. Chem., 2011)
    0.0166, %1/min TO BE FITTED
    0.0023, %nM/min (translation rate) TO BE FITTED
    973, %ratio concentration/fluorescence TO BE FITTED
    9, %autofluorescence;
    100, %nM/min (rate of enzymatic degradation) (Hersch, PNAS, 2004)
    75, %nM (dissociation constant of enzymatic degradation) (Hersch, PNAS, 2004)
    ]';

%specify the parameter bounds.
CanVary = [false,false,false,true,false,false,true,true,true,true,true,false,true,false,false, true, true, true, false,false,false];
lb = repmat(0,1,length(theta));
lb(~CanVary) = theta(~CanVary);
ub = repmat(Inf,1,length(theta));
ub(~CanVary) = theta(~CanVary);

%% Specify an initial guess for the state
%order of states in initial state vector is [s,m,s:m,c,p,z]
x0 = [1,1,1,1,1,10];
lb = [repmat(0,1,length(x0)), lb];
ub = [repmat(Inf,1,length(x0)), ub];
InitialGuess = [x0,theta];

%% plug into the wrapper function for the fitter
output = ParameterFit(InitialGuess,lb,ub);
