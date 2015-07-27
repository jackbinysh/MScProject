clc; clear variables;

%% Put in an initial guess, set upper and lower bounds, and get an answer out

% create a table containing bounds, and best guesses, for parameter values
RowNames = {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'beta';'c'};
VariableNames = {'LowerBounds','BestGuess','UpperBounds'};
LowerBounds = [10;10;10;0.001;0.1;0.01;0.01;0.0001;0.0001;10];
BestGuess = [98.0521504318544;560429.992864335;29126800.1065540;77.1270546961427;50.8028377485244;0.0817539315524263;61.7708875939270;0.00242232706183576;0.0217868661308902;421.538749015539];
UpperBounds = [1e3;1e6;1e8;100;100;100;100;10;1;10000];
Initialtheta = table(LowerBounds,BestGuess,UpperBounds,'RowNames',RowNames,'VariableNames',VariableNames);
clear RowNames VariableNames LowerBounds UpperBounds BestGuess

% create a table containing bounds, and best guesses, for initial state
RowNames = {'s0';'m0';'s_m0';'c0';'p0';'z0'};
VariableNames = {'LowerBounds','BestGuess','UpperBounds'};
LowerBounds = [0.1;0.1;0.1;0.1;0.1;9];
BestGuess = [117.4994; 987.8673;234.1104;751.4283;0.2720;9.5530];
UpperBounds = [1e3;1e3;1e3;1e3;1e3;11];
Initialx0 = table(LowerBounds,BestGuess,UpperBounds,'RowNames',RowNames,'VariableNames',VariableNames);
clear RowNames VariableNames LowerBounds UpperBounds BestGuess

% plug into the wrapper function for the fitter
 [xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta,Initialx0);

%% Plots etc
% the experimental data
Times = csvread('../data/InitialExperimentalData/time16_8.csv');
Data = csvread('../data/InitialExperimentalData/data16_8.csv');
% the solution found by optimisation
load('variablescmaes.mat','bestever');
% remember we have scaled this by our initial guess vector. So we scale
% back
bestever.x = bestever.x .* [Initialx0.BestGuess;Initialtheta.BestGuess];

% make tables out of the answer
RowNames = {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'beta';'c'};
VariableNames = {'BestGuess'};
theta = table(bestever.x(7:end),'RowNames',RowNames,'VariableNames',VariableNames);
RowNames ={'s0';'m0';'s_m0';'c0';'p0';'z0'};
x0 = table(bestever.x(1:6),'RowNames',RowNames,'VariableNames',VariableNames);

%%% REMEMBER TO COMMENT THIS OUT IF YOU DONT WANT IT
% if we want to use the initial guess for comparison of theta
%theta = Initialtheta(:,'BestGuess'); display('USING ORIGINAL THETA')

% now get the simulated predicition
ModelWithParams = @(t,x) RibodynamicsModel(t,x, theta);
[T,Prediction] = ode23s(ModelWithParams,Times,x0.BestGuess);

%plot the experimental data and simulation
plot(Times,Data); hold on;
plot(Times,Prediction(:,6),'linewidth', 3.5);



