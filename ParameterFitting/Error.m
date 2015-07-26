function Error = Error(p,Scale)

% scaling the parameter vector up
p = p .* Scale;
x0 = p(1:6);
theta = p(7:end);

% %reading in the experimental time series
 persistent Times ;
 if isempty(Times)
    Times = csvread('../data/InitialExperimentalData/time16_8.csv');
 end
 persistent Data;
 if isempty(Data)
    Data = csvread('../data/InitialExperimentalData/data16_8.csv');
 end

% % MATLAB
% % creating a function handle with the parameters built in, which
% % we then call in the ode solver
% %cd('./MatlabRibodynamics')
% %ModelWithParams = @(t,x) RibodynamicsModel(t,x, theta);
% % get the model prediction
% %[T,Prediction] = ode23s(ModelWithParams,Times,x0.BestGuess);
% %cd('..')
% 
% %MEX
cd('./MexRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
cd('..')

% compute the difference between them
Errors = Data - repmat(Prediction(:,6),1,length(Data));
Error = sum(sum(Errors.^2,1)); 
end
