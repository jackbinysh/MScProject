clc; clear variables;
format long g;

%Some bits and pieces of code for plotting
%This file contains various scraps of code I've used for plotting ,reading
%data in, etc.

% a workflow for reading in and analysing a single dataset


% the experimental data
Dataset = '13_9';
Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));
% remember we have scaled this by our initial guess vector. So we scale
% back

theta = bestever.x .* varargin{1}.Scale;

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
    name = strcat('./data/14_7_13_9_more/variablescmaes13_9_14_7_',num2str(i),'.mat');
    if exist( name,'file') == 2
        load(name,'bestever','varargin');
        f(end+1) = bestever.f;
        alltheta(:,end+1)  = bestever.x .* varargin{1}.Scale;
    end
end

filteredtheta = alltheta;
filteredf = f;

filteredtheta = alltheta(:, f < 5.55 & f~=0 );
filteredf = f(:, f < 5.55 & f~=0 );

filteredtheta = alltheta(:, f==min(f) );
filteredf = f(:, f==min(f) );

for i = 1:size(filteredtheta,2)
    theta = filteredtheta(:,i);
    x0 = InitialState(theta,0,1); 
    options = odeset('Jacobian',@Jacobian,'MaxStep',10); 
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'13_9');
    hold on;
    plot(Times,Data);
    predline = plot(Times,Prediction(:,2),'linewidth', 3.5);;
    expline = plot(Times,mean(Data,2),'--k','linewidth', 3.5);
end
% for 13_9
forcing = arrayfun(@(x) atc_input(x,'13_9'),Times)
SecondPartTime = 7*(60+60);
TimesIndex = Times<=SecondPartTime;
forcinglinefirstpart = plot(Times(TimesIndex),(forcing(TimesIndex) +14)/2, '--r','linewidth', 3.5); 
forcinglinesecondpart = plot(Times(~TimesIndex),(forcing(~TimesIndex) +14)/2, '--b','linewidth', 3.5); 
legend([predline,expline,forcinglinefirstpart,forcinglinesecondpart],'Prediction','Experimental mean','Forcing (60/60)', 'Forcing (60/30)')
% for 14_7
forcing = arrayfun(@(x) atc_input(x,'14_7'),Times)
forcingline = plot(Times,(forcing +14)/2, '--r','linewidth', 3.5); 
legend([predline,expline,forcingline],'Prediction','Experimental mean','Forcing')

for i = 1:length(filteredtheta)
    filteredx0(:,i) = InitialState(filteredtheta(:,i),1,1); 
end

%histograms of theta values 
figure
for i = 1:9
    subplot(3,3,i)
    histogram(filteredtheta(i,:),40)
end

% get the error of the whole timeseries


%histograms of function values 
figure
histogram(f,40)

% correlation between values
plot(filteredtheta(1,:),filteredtheta(2,:),'x')

%get the covraiance matrix
correlationmatrix = corrcoef(filteredtheta');
% remove the diagonals
correlationmatrix= correlationmatrix - diag(diag(correlationmatrix));
heatmap = HeatMap(correlationmatrix,'RowLabels',{'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'theta'},...
   'ColumnLabels',{'f_srna';'k_on';'k_off';'k_hyb';'d_m';'d_s';'mu';'beta';'theta'} )



%seeing what the parameters do
   
    theta(3) = 1000 %filteredtheta(3,1) - 0.01*filteredtheta(3,1);
    x0 = InitialState(theta,0,1);
    options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
    figure(1)
    hold on
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'13_9');
    plot(Times,Prediction(:,2),'linewidth', 3.5);
    figure(2)
    hold on
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'14_7');
    plot(Times,Prediction(:,6),'linewidth', 3.5);

% initiliase
    theta = alltheta(:,1);
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
theta = filteredtheta(:,2);
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
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
for i = 1:9
    hold on
    plot([0:1:max(Times)],Shat(1:10:end,i),markers{i}); xlim([120,240])
end
legend({'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'}')

figure
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
for i = 1:9
    hold on
    plot([0:1:max(Times)],Shat(1:10:end,i),markers{i}); xlim([120,240])
end
legend({'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'}')


% compute the co - linearity index
for i = 1:size(S,2)
    Z(:,i) = Shat(:,i)/norm(Shat(:,i));
end

%plot them normalised
figure
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
for i = 1:9
    hold on
    plot([0:1:max(Times)],abs(Z(1:10:end,i)),markers{i},'markers',5); xlim([120,240])
end
legend({'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'}')

%% analysis for argument about time scale forcing.

% get all the fixed points 

fixedpoints = [];
for i = 1:size(filteredtheta,2)
    theta = filteredtheta(:,i);
    x0 = InitialState(theta,1,1); 
    fixedpoints(:,end+1)  = x0;
end

lowfixedpoints = [];
for i = 1:size(filteredtheta,2)
    theta = filteredtheta(:,i);
    x0 = InitialState(theta,0,1); 
    lowfixedpoints(:,end+1)  = x0;
end

    options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
    figure
    hold on
    [T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,'13_9');
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    
% plot all state predictions from model, nomrlised
theta = filteredtheta(:,2);
x0 = InitialState(theta,0,1); 
options = odeset('Jacobian',@Jacobian,'MaxStep',10); 
[T,Prediction] = ode15s(@RibodynamicsModel,[0, max(Times)],x0,options,theta,'13_9');
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};
for i = 1:size(Prediction,2)
    hold on;
    % scale the prediction
    prediction = Prediction(:,i) - min( Prediction(:,i));
    prediction = prediction / max(prediction);
    predline = plot(T, prediction ,markers{i},'markers',9);
end
forcing = arrayfun(@(x) atc_input(x,'13_9'),T);
forcingline = plot(T,(forcing)/2, '--r','linewidth', 2); 
legend({'s';'m';'s:m';'c';'p';'z';'forcing';})
xlim([500,700]);

% analysis of fixed points etc

x = alltheta(8,:).*fixedpoints(2,:)+ alltheta(1,:).* alltheta(8,:).*fixedpoints(4,:);
filter = f >6 & f~=0;
y = alltheta(8,filter).*fixedpoints(2,filter)+ alltheta(1,filter).* alltheta(8,filter).*fixedpoints(4,filter);
hist(x,200)
plot(y,alltheta(7,filter),'x')
plot(f,x)
plot(alltheta(7,filter),fixedpoints(2,filter),'x')
corrcoef(alltheta(7,filter),fixedpoints(6,filter))

% plot mu against all state variable fixed points
for i = 1:9
    figure
    plot(alltheta(7,filter),fixedpoints(i,filter),'x')
    corrcoef(alltheta(7,filter),fixedpoints(i,filter))
end
% plot f against all state variable fixed points
for i = 1:9
    figure
    plot(f,fixedpoints(i,:),'x')
end


%%% making the heatmaps

f=[];alltheta=[];
for i = 1:200
    %name = strcat('./data/14_7_more/variablescmaes14_7_',num2str(i),'.mat');
    name = strcat('./data/13_9_more/variablescmaes13_9_',num2str(i),'.mat');
    if exist( name,'file') == 2
        load(name,'bestever','varargin');
        f(end+1) = bestever.f;
        alltheta(:,end+1)  = bestever.x .* varargin{1}.Scale;
    end
end

filteredtheta = alltheta;
filteredf = f;

% correlation between values
plot(filteredtheta(1,:),filteredtheta(2,:),'x')

%get the covraiance matrix
correlationmatrix = corrcoef(filteredtheta');
% remove the diagonals
correlationmatrix= correlationmatrix - diag(diag(correlationmatrix));
%heatmapone = HeatMap(correlationmatrix,'DisplayRange', 0.8,'RowLabels',{'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'theta'},...
%   'ColumnLabels',{'f_srna';'k_on';'k_off';'k_hyb';'d_m';'d_s';'mu';'beta';'theta'} )
imagesc(correlationmatrix,[-0.8,0.8])
set(gca,'XTickLabel',{'f_{srna}';'k_{on}';'k_{off}';'k_{hyb}';'delta_{m}';'delta_{s}';'mu';'beta';'theta'})
set(gca,'YTickLabel',{'f_{srna}';'k_{on}';'k_{off}';'k_{hyb}';'delta_{m}';'delta_{s}';'mu';'beta';'theta'})

f=[];alltheta=[];
for i = 1:200
    name = strcat('./data/14_7_more/variablescmaes14_7_',num2str(i),'.mat');
    if exist( name,'file') == 2
        load(name,'bestever','varargin');
        f(end+1) = bestever.f;
        alltheta(:,end+1)  = bestever.x .* varargin{1}.Scale;
    end
end

filteredtheta = alltheta;
filteredf = f;

%get the covraiance matrix
correlationmatrix = corrcoef(filteredtheta');
% remove the diagonals
correlationmatrix= correlationmatrix - diag(diag(correlationmatrix));
%heatmaptwo = HeatMap(correlationmatrix,'DisplayRange', 0.8,'RowLabels',{'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'theta'},...
%   'ColumnLabels',{'f_srna';'k_on';'k_off';'k_hyb';'d_m';'d_s';'mu';'beta';'theta'} )
figure
imagesc(correlationmatrix,[-0.8,0.8])
set(gca,'XTickLabel',{'f_{srna}';'k_{on}';'k_{off}';'k_{hyb}';'delta_{m}';'delta_{s}';'mu';'beta';'theta'})
set(gca,'YTickLabel',{'f_{srna}';'k_{on}';'k_{off}';'k_{hyb}';'delta_{m}';'delta_{s}';'mu';'beta';'theta'})

for i = 1:9
    figure(1)
    h = subplot(3,3,i);
    h.FontSize = 20;
    h.XLabel.FontSize = 20
    h.YLabel.FontSize = 20
    Pos = h.Position;
    figure(2);
    h = subplot(3,3,i);
    h.Position = Pos;
    h.FontSize = 20;
    h.XLabel.FontSize = 20
    h.YLabel.FontSize = 20
end


