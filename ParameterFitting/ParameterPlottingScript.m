%% Some bits and pieces of code for plotting
% the experimental data
Times = csvread('./data/InitialExperimentalData/time16_8.csv');
Data = csvread('./data/InitialExperimentalData/data16_8.csv');
% the solution found by optimisation
load('variablescmaes.mat');
% remember we have scaled this by our initial guess vector. So we scale
% back

theta = out.solutions.recentbest.x .* varargin{1}.Scale;
x0 = InitialState(theta); 

%MEX
cd('./MexRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
cd('..')

%plot the experimental data and simulation
hold off
plot(Times,Data); hold on;
plot(Times,Prediction(:,6),'linewidth', 3.5);