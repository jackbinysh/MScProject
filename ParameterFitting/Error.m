function Error = Error(p)

% the first 6 parameters are the initial state
x0 = p(1:6); 
% the rest are the model parameters 
theta = p(7:end);

%reading in the experimental time series
persistent Times ;
Times = csvread('../data/InitialExperimentalData/time16_8.csv');
persistent Data;
Data = csvread('../data/InitialExperimentalData/data16_8.csv');

% craftily creating a function handle with the parameters built in, which
% we then call in the ode solver
ModelWithParams = @(t,x) RibodynamicsModel(t,x, theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
% compute the difference between them
Errors = Data - repmat(Prediction(:,6),1,length(Data));
Error = sum(sum(Errors.^2,1)); 
end
