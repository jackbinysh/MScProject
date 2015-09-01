function Error = Error(theta,Varargin)

% scaling the parameter vector up
Scale = Varargin.Scale;
theta = theta .* Scale;
Data = Varargin.Data;


%get our initial state from the models fixed point, given the parameters
% the fixed point is solved for with IPTG on.
x0 = InitialState(theta,0,1);
% Set the ODE solver options. Provide a jacobian, and ensure the solver doesn't
% step over a period of oscillation by setting a maximum step size.
% http://uk.mathworks.com/help/matlab/ref/odeset.html
options = odeset('Jacobian',@Jacobian,'MaxStep',10);
% now, for each dataset, get the prediction, find its error, add it to the
% running total

% note, since we are unsure of our initial state, but we reach the
% attractor quickly, I'll just discard the first 10% of the data.
Error = 0;
for i = 1:length(Data(1,:))
    Dataset = Data{1,i};
    Times = Data{2,i};
    Fluorescence = Data{3,i};
    
    Startpoint = ceil(length(Times)/10); % setting the amount to discard]
    SecondPartTime = 7* (60 + 60); % 7 initial periods
    BelowTime = Times < SecondPartTime;
    Endpoint = max(find(BelowTime));
 
    % run the simulation
    [~,Prediction] = ode15s(@RibodynamicsModel,Times,x0,options,theta,Dataset);
    
    % sometimes the integrator will fail due to step size tolerance. these
    % are bad points, and we just set the error to inf for them. They
    % occur in the initial guessing phase of the algorithm.
    try
        Error = Error + sum((Fluorescence(Startpoint:Endpoint) - Prediction(Startpoint:Endpoint,6)).^2);
    catch
        Error = Inf;
    end
        
end
end
