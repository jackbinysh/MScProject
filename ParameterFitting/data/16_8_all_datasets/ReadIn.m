% a bit of code to analyse all the different time series
for i = 1:145
    name = strcat('./data/16_8_all_datasets/variablescmaes',num2str(i),'.mat');
    load(name,'bestever','varargin');
    alltheta(:,i)  = bestever.x .* varargin{1}.Scale;
    theta = alltheta(:,i);
    x0 = InitialState(theta); 

    %MATLAB
    cd('./MatlabRibodynamics')
    ModelWithParams = @(t,x) RibodynamicsModel(t,x', theta);
    % get the model prediction
    [T,Prediction] = ode15s(ModelWithParams,Times,x0);
    cd('..')
    hold on;
    %plot the experimental data and simulation
    plot(Times,Prediction(:,6),'linewidth', 3.5);
    
    
end

theta = mean(alltheta,2);
