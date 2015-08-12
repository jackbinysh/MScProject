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
    
    % constructing the struct of additional arguments to pass to the objective
    Varargin = struct('Scale',Scale,'Dataset',Dataset);
    
    % setting some display etc. options
    opts.DispModulo = 1;
    opts.Restarts = 2;
    opts.SaveFilename = strcat('variablescmaes',Dataset,num2str(SaveNumber),'.mat');
    opts.LogFilenamePrefix = strcat('outcmaes',Dataset,num2str(SaveNumber));
    
    
    % passing everyting into the algorithm
    [xmin,fmin,counteval,stopflag,out] = cmaes('Error',Initialtheta,[],opts,Varargin);
end





