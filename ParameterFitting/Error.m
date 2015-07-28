function Error = Error(p,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
p = p .* Scale;
x0 = p(1:6);
theta = p(7:end);

%reading in the experimental time series
 persistent Times Data MeanFluorescence;
 if (isempty(Times))
    Dataset = Varargin.Dataset;
    Times = csvread(strcat('../data/InitialExperimentalData/time',Dataset,'.csv'));
    Data = csvread(strcat('../data/InitialExperimentalData/data',Dataset,'.csv'));
    MeanFluorescence = mean(Data,2);
 end

%MEX
cd('./MexRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
% get the model prediction
[T,Prediction] = ode23s(ModelWithParams,Times,x0);
cd('..')

Error = sum((MeanFluorescence - Prediction(:,6)).^2);
end
