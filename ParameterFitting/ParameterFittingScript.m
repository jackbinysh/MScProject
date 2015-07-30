clc; clear variables;
format long g;

% decide which dataset we are using
Dataset = '16_8';

% decide which run we use: if we leave it empty, an average is taken
RunNumber = 1;

% Initial Guess and bounds
% list of variables {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'beta';'c'};
InitialthetaLB = [10;10;10;0.001;0.1;0.01;0.01;0.0001;0.0001;10];
Initialtheta = [98.0521504318544;560429.992864335;29126800.1065540;77.1270546961427;50.8028377485244;0.0817539315524263;61.7708875939270;0.00242232706183576;0.0217868661308902;421.538749015539];
InitialthetaUB = [1e3;1e6;1e8;100;100;100;100;10;1;10000];

% plug into the wrapper function for the fitter
[xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta, InitialthetaLB, InitialthetaUB, Dataset,RunNumber);