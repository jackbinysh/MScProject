function Error = Error(theta,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
theta = theta .* Scale;
Dataset = Varargin.Dataset;

%reading in the experimental time series
 persistent Times Data Fluorescence Startpoint;
 if (isempty(Times))
    Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
    Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));
    Fluorescence = mean(Data,2);
    % we are not sure of the initial state, but the model rapidly finds its
    % fixed point. We just discard the first 5% of the data
    Startpoint = ceil(length(Times/20));
 end

%get our initial state from the models fixed point, given the parameters
% the fixed point is solved for with IPTG on.
cd('./model1')
x0 = InitialState(theta,0,1);

% Set the ODE solver options. Provide a jacobian, and ensure the solver doesn't
% step over a period of oscillation by setting a maximum step size.
% http://uk.mathworks.com/help/matlab/ref/odeset.html
options = odeset('Jacobian',@Jacobian,'MaxStep',10); 
[T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,Dataset);
cd('..')

Error = sum((Fluorescence(Startpoint:end) - Prediction(StartPoint:end,6)).^2);
end
