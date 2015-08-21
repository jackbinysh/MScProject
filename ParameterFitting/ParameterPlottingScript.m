clc; clear variables;
format long g;

%Some bits and pieces of code for plotting
%This file contains various scraps of code I've used for plotting ,reading
%data in, etc.

% a workflow for reading in and analysing a single dataset


% the experimental data
Dataset = '14_7';
Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));
% remember we have scaled this by our initial guess vector. So we scale
% back

theta = bestever.x .* varargin{1}.Scale;

cd('./model1')
x0 = InitialState(theta,0,1); 
options = odeset('Jacobian',@Jacobian,'MaxStep',10); 
[T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,Dataset);
cd('..')

hold on;
plot(Times,Data);
plot(Times,Prediction(:,6),'linewidth', 3.5);
plot(Times,mean(Data,2),'--k','linewidth', 3.5);
hold off

% for many curves

% the experimental data
Dataset = '13_9';
Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));

f=[];alltheta=[];
for i = 1:100
    name = strcat('./data/interim/variablescmaes13_9_14_7_',num2str(i),'.mat');
    if exist( name,'file') == 2
        load(name,'bestever','varargin');
        f(end+1) = bestever.f;
        alltheta(:,end+1)  = bestever.x .* varargin{1}.Scale;
    end
end

filteredtheta = alltheta(:, f < 30 & f~=0 );
filteredf = f(:, f < 30 & f~=0 );

for i = 1:size(filteredtheta,2)
    theta = filteredtheta(:,i);
    x0 = InitialState(theta,0,1); 
    options = odeset('Jacobian',@Jacobian,'MaxStep',10); 
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'13_9');
    hold on;
    plot(Times,Data);
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    plot(Times,mean(Data,2),'--k','linewidth', 3.5);
end


for i = 1:length(filteredtheta)
    filteredx0(:,i) = InitialState(filteredtheta(:,i),1,1); 
end

%histograms of theta values 
figure
for i = 1:9
    subplot(3,3,i)
    histogram(filteredtheta(i,:),20)
end

% correlation between values
plot(filteredtheta(1,:),filteredtheta(2,:),'x')

%seeing what the parameters do
   
    theta(7) = filteredtheta(7,1) - 0.01*filteredtheta(7,1);
    x0 = InitialState(theta,0,1);
    options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
    figure(1)
    hold on
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'13_9');
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    figure(2)
    hold on
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'14_7');
    plot(Times,Prediction(:,6),'linewidth', 3.5);

% initiliase
    theta = filteredtheta(:,1);
    x0 = InitialState(theta,0,1); 
    options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
    figure(1)
    hold on
    Dataset = '13_9';
    Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
    Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'13_9');
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    plot(Times,mean(Data,2),'--k','linewidth', 3.5);
    figure(2)
    hold on
    Dataset = '14_7';
    Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
    Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'14_7');
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    plot(Times,mean(Data,2),'--k','linewidth', 3.5);


% a little code snippet to plot the forcing
forcing = arrayfun(@(x) atc_input(x,'13_9'),Times)

for i = 1:9
    subplot(3,3,i);
    histogram(filteredtheta(i,:),20);
end


figure
for i = 1:9
    subplot(3,3,i)
    plot([0:0.1:max(Times)],S(i,:)); xlim([120,240])
end

figure
for i = 1:9
    hold on;
    plot([0:0.1:max(Times)],x(i,:)/max(abs(x(i,:)))); xlim([120,240])
end
legend({'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'}')

%% A sensitvity Analysis

% get the sensitivity matrix
theta = filteredtheta(:,1);
S = SensitivityMatrix(theta,'13_9');

% scale them by the guess
for i = 1:size(S,2)
    Shat(:,i) = S(:,i)*theta(i);
end

% plot them
figure
for i = 1:9
    subplot(3,3,i)
    plot([0:0.1:max(Times)],Shat(:,i)); xlim([120,240])
end


% plot them all on the same graph
figure
for i = 1:9
    hold on
    plot([0:0.1:max(Times)],Shat(:,i)); xlim([120,240])
end

%plot them normalised
figure
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
for i = 1:9
    hold on
    plot([0:1:max(Times)],Shat(1:10:end,i),markers{i}); xlim([120,240])
end
legend({'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'}')


% compute the co - linearity index
for i = 1:size(S,2)
    Z(:,i) = S(:,i)/norm(S(:,i));
end

%plot them normalised
figure
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
for i = 1:9
    hold on
    plot([0:1:max(Times)],abs(Z(1:10:end,i)),markers{i}); xlim([120,240])
end
legend({'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'}')

    
