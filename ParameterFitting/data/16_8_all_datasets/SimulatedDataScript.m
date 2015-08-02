% a script to generate some artificial data, which 
% we can test our algorithm on

% artificial parameters
theta = [100;560000;29126800;77;50;0.08;62;0.003;0.003;422];
x0 = InitialState(theta);

% generate model prediction
Times = csvread('./data/InitialExperimentalData/time16_8.csv');

%MEX
cd('./MexRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
cd('..')

Prediction = Prediction(:,6);
% now fuzz the time series with gaussian noise.
%First set a scale for the noise
Scale = 0.2*(max(Prediction) - mean(Prediction));

% now generate noisy paths
NumSeries = 145;
Noise = normrnd(0,repmat(Scale,length(Prediction),NumSeries));
SimulatedData = repmat(Prediction,1,NumSeries) + Noise;

% create a table containing bounds, and best guesses, for starting parameter values
RowNames = {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'beta';'c'};
VariableNames = {'LowerBounds','BestGuess','UpperBounds'};
LowerBounds = [10;10;10;0.001;0.1;0.01;0.01;0.0001;0.0001;10];
UpperBounds = [1e3;1e6;1e8;100;100;100;100;10;1;10000];
% we get the best guess from uniform guesses across the plausible range
BestGuess = unifrnd(LowerBounds,UpperBounds);
Initialtheta = table(LowerBounds,BestGuess,UpperBounds,'RowNames',RowNames,'VariableNames',VariableNames);
clear RowNames VariableNames LowerBounds UpperBounds BestGuess

% create a table containing bounds, and best guesses, for initial state
RowNames = {'s0';'m0';'s_m0';'c0';'p0';'z0'};
VariableNames = {'LowerBounds','BestGuess','UpperBounds'};
LowerBounds = [0.1;0.1;0.1;0.1;0.1;9];
UpperBounds = [1e3;1e3;1e3;1e3;1e3;11];
BestGuess = unifrnd(LowerBounds,UpperBounds);
Initialx0 = table(LowerBounds,BestGuess,UpperBounds,'RowNames',RowNames,'VariableNames',VariableNames);
clear RowNames VariableNames LowerBounds UpperBounds BestGuess

% plug into the wrapper function for the fitter
 [xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta,Initialx0);

