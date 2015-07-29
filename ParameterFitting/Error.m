function Error = Error(theta,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
theta = theta .* Scale;

%get our initial state from the models fixed point, given the parameters
x0 = InitialState(theta);

%reading in the experimental time series
 persistent Times Data MeanFluorescence;
 if (isempty(Times))
    Dataset = Varargin.Dataset;
    Times = csvread(strcat('./data/InitialExperimentalData/time',Dataset,'.csv'));
    Data = csvread(strcat('./data/InitialExperimentalData/data',Dataset,'.csv'));
    MeanFluorescence = mean(Data,2);
 end

%MEX
cd('./MatlabRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
% get the model prediction
[T,Prediction] = ode15s(ModelWithParams,Times,x0);
cd('..')

Error = sum((MeanFluorescence - Prediction(:,6)).^2);
end
