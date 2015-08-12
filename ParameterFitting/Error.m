function Error = Error(theta,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
theta = theta .* Scale;
Dataset = Varargin.Dataset;

%reading in the experimental time series
 persistent Times Data Fluorescence;
 if (isempty(Times))
    Times = csvread(strcat('./data/CleanedData/time',Dataset,'.csv'));
    Data = csvread(strcat('./data/CleanedData/data',Dataset,'.csv'));
    Fluorescence = mean(Data,2);
 end

%get our initial state from the models fixed point, given the parameters
% the fixed point is solved for with IPTG on.
cd('./model1')
x0 = InitialState(theta,0,1);

% Set the ODE solver options. Provide a jacobian, and ensure the solver doesn't
% step over a period of oscillation by setting a maximum step size.
% http://uk.mathworks.com/help/matlab/ref/odeset.html
options = odeset('Jacobian',@Jacobian,'MaxStep',30); 
[T,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,Dataset);
cd('..')

Error = sum((Fluorescence - Prediction(:,6)).^2);
end
