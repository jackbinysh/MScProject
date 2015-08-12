clc; clear variables;
format long g;

%Some bits and pieces of code for plotting
%This file contains various scraps of code I've used for plotting ,reading
%data in, etc.

% a workflow for reading in and analysing a single dataset

% the solution found by optimisation
load('variablescmaes.mat');

% the experimental data
Dataset = '16_8';
Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));
% remember we have scaled this by our initial guess vector. So we scale
% back

theta = bestever.x .* varargin{1}.Scale;

% plotting the original model with its original parameters
% order should be {'f_srna';'k_on';'k_off';'k_hyb';'k_matur';'delta_m';'delta_s';'mu';'beta';'c'};
theta = [435,1e5,1e7,0.01,0.175,0.1,0.0166,0.0023,973]'

x0 = InitialState(theta,0,1); 
options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
[T,Prediction] = ode15s(@RibodynamicsModel,[0,Times(end)],x0,options,theta,Dataset);


%plot the experimental data and simulation
hold on;
plot(Times,Data);
plot(T,Prediction(:,6),'linewidth', 3.5);
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

for i = 1:200
    name = strcat('./data/16_8_mean_new_model_manymoretimes/variablescmaes',num2str(i),'.mat');
    load(name,'bestever','varargin');
    f(i) = bestever.f;
    alltheta(:,i)  = bestever.x .* varargin{1}.Scale;
end

filteredtheta = alltheta(:,f<3);

for i = 1:length(filteredtheta)
    theta = filteredtheta(:,i);
    x0 = InitialState(theta,0,1); 
    options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'16_8');
    hold on;
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    onstate = InitialState(theta,1,1);
    plot(Times,repmat(onstate(6),1,length(Times)));
    plot(Times,Data);

end

for i = 1:length(filteredtheta)
    filteredx0(:,i) = InitialState(filteredtheta(:,i),1,1); 
end

%histograms of theta values 
figure
for i = 1:9
    subplot(3,3,i)
    histogram(filteredtheta(i,:)/mean(filteredtheta(i,:)),20)
end

% correlation between values
plot(filteredtheta(1,:),filteredtheta(2,:),'x')

%seeing what the parameters do
    theta = theta + [100 0 0 0 0 0 0 0 0]'
    x0 = InitialState(theta,0,1);
    onstate = InitialState(theta,1,1);
    options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'16_8');
    hold on
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    plot(Times,repmat(onstate(6),1,length(Times)));

% reset
    theta = filteredtheta(:,1);
    x0 = InitialState(theta,0,1); 
    onstate = InitialState(theta,1,1);
    options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'16_8');
    hold on
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    plot(Times,repmat(onstate(6),1,length(Times)));


% a little code snippet to plot the forcing
forcing = arrayfun(@(x) atc_input(x,'13_9'),Times)

for i = 3:9
    histogram(filteredtheta(i,:),500)
    export_fig -pdf feval(num2str(i))
end

    
