%% Some bits and pieces of code for plotting
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