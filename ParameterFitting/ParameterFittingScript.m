clc; clear variables;

% create a table containing bounds, and best guesses, for parameter values
RowNames = {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'Beta';'c'};
VariableNames = {'LowerBounds','BestGuess','UpperBounds'};
LowerBounds = [0;0;0;0;0;0;0;0;0;0];;
BestGuess = [435;1e5; 1e7; 0.01; 0.175; 0.1; 0.1; 0.0166; 0.0023; 973]; 
UpperBounds = [Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf];
theta = table(LowerBounds,BestGuess,UpperBounds,'RowNames',RowNames,'VariableNames',VariableNames);
clear RowNames VariableNames LowerBounds UpperBounds BestGuess

% create a table containing bounds, and best guesses, for initial state
RowNames = {'s0';'m0';'s_m0';'c0';'p0';'z0'};
VariableNames = {'LowerBounds','BestGuess','UpperBounds'};
LowerBounds = [0;0;0;0;0;8];
BestGuess = [10;10;10;10;10;9]; 
UpperBounds = [Inf;Inf;Inf;Inf;Inf;11];
x0 = table(LowerBounds,BestGuess,UpperBounds,'RowNames',RowNames,'VariableNames',VariableNames);
clear RowNames VariableNames LowerBounds UpperBounds BestGuess

%% plug into the wrapper function for the fitter
[History,output] = ParameterFit(theta,x0);
clear theta, x0;

%% plot the prediction with these variables
x0  = output.x(end,1:6);
theta = struct('f_srna',output.x(7),'k_on',output.x(8),'k_off', output.x(9),'k_hyb', output.x(10),'delta_m', output.x(11),'delta_s',...
output.x(12),'delta_c',output.x(13),'mu',output.x(14),'Beta',output.x(15),'c',output.x(16));
% we then call in the ode solver
ModelWithParams = @(t,x) RibodynamicsModel(t,x, theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);


