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
% order should be {'mu';'F';'c'};
InitialthetaLB = [0.001; 10 ; 100];
InitialthetaUB = [0.1; 1000 ; 10000];

% for an initial theta, we uniformly select a point in this hypercube
Initialtheta = InitialthetaLB + (InitialthetaUB-InitialthetaLB ).*rand(length(InitialthetaLB),1) % i want it displayed

% plug into the wrapper function for the fitter
[xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta, InitialthetaLB, InitialthetaUB, Dataset,SaveNumber);