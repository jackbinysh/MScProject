clc; clearvars -except SaveNumber Dataset;
format long g;

% decide which dataset we are using,and 
% Give a save number, if we want to do multiple runs.
% im often setting this outside of the script, so I'll only set this now if
% its empty
if(~exist('SaveNumber','var')) SaveNumber = [1]; end
if(~exist('Dataset','var')) Dataset = {'13_9'}; end;

% put a random seed in to ensure different cluster runs don't start at the
% same point 
c = clock; secs = c(6);
seed = SaveNumber + secs + cputime % display it to check its always different
rng(seed, 'twister');

% Initial Guess and bounds. We ballpark these for now, based on
% the histograms of values we have gotten from earlier analyses.
% order should be {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'};
InitialthetaLB = [0.1;0.1;0.1;0.1;1;0.1;0.001;0.0001;300]
InitialthetaUB = [1e4;1e7;1e8;10000;10000;1000;5;100;2000]

% for an initial theta, we uniformly select a point in this hypercube
%Initialtheta = InitialthetaLB + (InitialthetaUB-InitialthetaLB ).*rand(length(InitialthetaLB),1) % i want it displayed

Initialtheta = [1476.09315060305;91027.8669948397;68425985.7997264;1;451.715025075103;1;0.0500000000000000;10;507.453643463837]

% plug into the wrapper function for the fitter
[xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta, InitialthetaLB, InitialthetaUB, Dataset,SaveNumber);