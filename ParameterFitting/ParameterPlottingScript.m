clc; clear variables;
format long g;

%Some bits and pieces of code for plotting
%This file contains various scraps of code I've used for plotting ,reading
%data in, etc.

% a workflow for reading in and analysing a single dataset

% the solution found by optimisation
load('variablescmaes.mat');

% the experimental data
Dataset = '13_9';
Times = csvread(strcat('./data/InitialExperimentalData/time',Dataset,'.csv'));
Data = csvread(strcat('./data/InitialExperimentalData/data',Dataset,'.csv'));
% remember we have scaled this by our initial guess vector. So we scale
% back

theta = bestever.x .* varargin{1}.Scale;
x0 = InitialState(theta,0,0); 

%MATLAB
cd('./MatlabRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta, Dataset);
% get the model prediction
options = odeset('MaxStep',1);
sol = ode15s(ModelWithParams,Times,x0,options);
[T,Prediction ]= ode15s(ModelWithParams,Times,x0,options);
cd('..')

%plot the experimental data and simulation
hold on;
plot(Times,Data);
plot(Times,Prediction(:,6),'linewidth', 3.5);
hold off

% the residuals
figure(2);
Residuals = Data - repmat(Prediction(:,6),1,length(Data));
hold off
plot(T,Residuals);

% a bit of code to analyse all the different time series
clear all;

Times = csvread('./data/InitialExperimentalData/time16_8.csv');
Data = csvread('./data/InitialExperimentalData/data16_8.csv');

for i = 1:145
    name = strcat('./data/16_8_all_datasets/variablescmaes',num2str(i),'.mat');
    load(name,'bestever','varargin');
    f(i) = bestever.f;
    alltheta(:,i)  = bestever.x .* varargin{1}.Scale;
end

for i = 1:length(alltheta)
    theta = alltheta(:,i);
    x0 = InitialState(theta); 

    %MATLAB
    cd('./MatlabRibodynamics')
    ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
    % get the model prediction
    [T,Prediction] = ode15s(ModelWithParams,Times,x0);
    cd('..')
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    plot(Times,Data);
    hold on;
end

% a little code snippet to plot the forcing
forcing = arrayfun(@(x) atc_input(x,'13_9'),Times)

    
