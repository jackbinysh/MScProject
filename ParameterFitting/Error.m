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
%MEX
cd('./MexRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
cd('..')

% compute the difference between them
Errors = Data - repmat(Prediction(:,6),1,length(Data));
Error = sum(sum(Errors.^2,1)); 
end
