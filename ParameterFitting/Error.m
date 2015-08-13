function Error = Error(theta,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
theta = theta .* Scale;
Data = Varargin.Data;

% choose the model we want to use
cd('./model1')
%get our initial state from the models fixed point, given the parameters
% the fixed point is solved for with IPTG on.
x0 = InitialState(theta,0,1);
% Set the ODE solver options. Provide a jacobian, and ensure the solver doesn't
% step over a period of oscillation by setting a maximum step size.
% http://uk.mathworks.com/help/matlab/ref/odeset.html
options = odeset('Jacobian',@Jacobian,'MaxStep',10);
% now, for each dataset, get the prediction, find its error, add it to the
% running total
Error = 0;
for i = 1:length(Data(1,:))
    Dataset = Data{1,i};
    Times = Data{2,i};
    Fluorescence = Data{3,i};
    
    [~,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,Dataset);
    Error = Error + sum((Fluorescence - Prediction(:,6)).^2);
end
cd('..')
end
