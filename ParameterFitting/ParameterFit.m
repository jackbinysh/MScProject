function [xmin,fmin,counteval,stopflag,out] = ParameterFit(Initialtheta,InitialthetaLB,InitialthetaUB,Dataset,SaveNumber)

    % reading off the initial guess and the bounds
    opts.LBounds = InitialthetaLB;  
    opts.UBounds = InitialthetaUB;
    
    % The cmaes wants all parameters to be on the same scale. I take our
    % initial Guess as setting the scale for the problem
    Scale =  Initialtheta;
    opts.LBounds = opts.LBounds./Scale;
    opts.UBounds = opts.UBounds./Scale;
    Initialtheta = Initialtheta./Scale;
    
    % loading in the datasets we want
    for i = 1:length(Dataset)
        Data{1,i} = Dataset{i};
        Data{2,i} = csvread(strcat('./data/CleanedData/time',Dataset{i},'.csv'));
        Data{3,i} = mean(csvread(strcat('./data/CleanedData/data',Dataset{i},'.csv')),2);
    end
    
    % constructing the struct of additional arguments to pass to the objective
    Varargin = struct('Scale',Scale,'Data',{Data});
    
    % setting some display etc. options
    opts.DispModulo = 1;
    opts.Restarts = 2;
    opts.SaveFilename = strcat('variablescmaes',strjoin(Dataset,'_'),'_',num2str(SaveNumber),'.mat');
    opts.LogFilenamePrefix = strcat('outcmaes',strjoin(Dataset,'_'),'_',num2str(SaveNumber));
    
    
    % passing everyting into the algorithm
    [xmin,fmin,counteval,stopflag,out] = cmaes('Error',Initialtheta,[],opts,Varargin);
end





