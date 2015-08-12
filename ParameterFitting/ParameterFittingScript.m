clc; clear variables;
format long g;

% decide which dataset we are using
Dataset = '16_8';

% decide which run we use: if we leave it empty, an average is taken
RunNumber = [];

% Initial Guess and bounds. Thesse  are currently set to the best guess
% from GArun_30_07_2015.
%  % order should be {'f_srna';'k_on';'k_off';'k_hyb';'k_matur';'delta_m';'delta_s';'mu';'beta';'c'};
InitialthetaLB = [1;1000;100000;1;1;1;0.00001;0.0001;1];
InitialthetaUB = [1e3;1e6;1e8;1000;1000;1000;1;10;10000];
Initialtheta = [500;5e5;5e7;500;500;500;0.02;1;450];

% plug into the wrapper function for the fitter
[xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta, InitialthetaLB, InitialthetaUB, Dataset,RunNumber);