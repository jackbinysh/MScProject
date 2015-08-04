function Error = Error(theta,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
theta = theta .* Scale;

%get our initial state from the models fixed point, given the parameters
x0 = InitialState(theta,0,0);

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

%MATLAB
cd('./MatlabRibodynamics')
ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta,Varargin.Dataset);
% Set the ODE solver options. Provide a jacobian, and ensure the solver doesn't
% step over a period of oscillation by setting a maximum step size.
% http://uk.mathworks.com/help/matlab/ref/odeset.html

options = odeset('Jacobian',@jacobian,'MaxStep',30); 
[T,Prediction] = ode15s(ModelWithParams,Times,x0,options);
cd('..')

Error = sum((Fluorescence - Prediction(:,6)).^2);
end
