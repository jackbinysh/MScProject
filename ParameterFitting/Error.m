function Error = Error(theta,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
theta = theta .* Scale;

%get our initial state from the models fixed point, given the parameters
x0 = InitialState(theta);

%reading in the experimental time series
 persistent Times Data Fluorescence;
 if (isempty(Times))
    Times = csvread(strcat('./data/InitialExperimentalData/time',Varargin.Dataset,'.csv'));
    Data = csvread(strcat('./data/InitialExperimentalData/data',Varargin.Dataset,'.csv'));
    
    %if we've specified a specific run we are interested in, thats our data
    if(~isempty(Varargin.RunNumber))
        Fluorescence = Data(:,Varargin.RunNumber);
    else
    % if its empty, take the mean of all the data
        Fluorescence = mean(Data,2);
    end
 end

%MEX
cd('./MatlabRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
% get the model prediction
[T,Prediction] = ode15s(ModelWithParams,Times,x0);
cd('..')

Error = sum((Fluorescence - Prediction(:,6)).^2);
end
