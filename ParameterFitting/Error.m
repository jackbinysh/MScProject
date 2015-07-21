function Error = Error(p)

% the first 6 parameters are a guess for the initial state
x0 = p(1:6); 
% the rest are a guess for the model parameters 
theta = struct('f_srna',p(7),'k_on',p(8),'k_off', p(9),'k_hyb', p(10),'delta_m', p(11),'delta_s',...
p(12),'delta_c',p(13),'mu',p(14),'Beta',p(15),'c',p(16));

%reading in the experimental time series
persistent Times ;
Times = csvread('../data/InitialExperimentalData/time16_8.csv');
persistent Data;
Data = csvread('../data/InitialExperimentalData/data16_8.csv');

% creating a function handle with the parameters built in, which
% we then call in the ode solver
ModelWithParams = @(t,x) RibodynamicsModel(t,x, theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
% compute the difference between them
Errors = Data - repmat(Prediction(:,6),1,length(Data));
Error = sum(sum(Errors.^2,1)); 
end
