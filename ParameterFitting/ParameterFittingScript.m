clc; clear variables;

%% decide which dataset we are using
Dataset = '16_8';

%% create 2 tables, giving bounds on the initial parameters and state
% create a table containing bounds, and best guesses, for parameter values
RowNames = {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'beta';'c'};
VariableNames = {'LowerBounds','BestGuess','UpperBounds'};
LowerBounds = [10;10;10;0.001;0.1;0.01;0.01;0.0001;0.0001;10];
BestGuess = [98.0521504318544;560429.992864335;29126800.1065540;77.1270546961427;50.8028377485244;0.0817539315524263;61.7708875939270;0.00242232706183576;0.0217868661308902;421.538749015539];
UpperBounds = [1e3;1e6;1e8;100;100;100;100;10;1;10000];
Initialtheta = table(LowerBounds,BestGuess,UpperBounds,'RowNames',RowNames,'VariableNames',VariableNames);
clear RowNames VariableNames LowerBounds UpperBounds BestGuess

%% plug into the wrapper function for the fitter
 [xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta,Dataset);