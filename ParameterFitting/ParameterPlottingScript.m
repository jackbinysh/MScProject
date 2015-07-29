%% Some bits and pieces of code for plotting
% the experimental data
Times = csvread('../data/InitialExperimentalData/time16_8.csv');
Data = csvread('../data/InitialExperimentalData/data16_8.csv');
% the solution found by optimisation
load('variablescmaes.mat','bestever');
% remember we have scaled this by our initial guess vector. So we scale
% back
BestGuess = [98.0521504318544;560429.992864335;29126800.1065540;77.1270546961427;50.8028377485244;0.0817539315524263;61.7708875939270;0.00242232706183576;0.0217868661308902;421.538749015539];
theta = bestever.x .*BestGuess;
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