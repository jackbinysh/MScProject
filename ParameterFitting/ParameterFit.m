function [xmin,fmin,counteval,stopflag,out] = ParameterFit(theta,Dataset)

    % reading off the initial guess and the bounds
    BestGuess = theta.BestGuess;
    opts.LBounds = theta.LowerBounds;  
    opts.UBounds = theta.UpperBounds;
    
    % The cmaes wants all parameters to be on the same scale. I take our
    % initial Guess as setting the scale for the problem
    Scale =  BestGuess;
    opts.LBounds = opts.LBounds./Scale;
    opts.UBounds = opts.UBounds./Scale;
    BestGuess = BestGuess./Scale;
    
    % setting an initial spread of points, using 1/3rd of the search region as recommended in cmaes script
    Sigma = (opts.UBounds - opts.LBounds)/3; 
    
    % constructing the struct of additional arguments to pass to the objective
    Varargin = struct('Scale',Scale,'Dataset',Dataset);
    
    % setting some display etc. options
    opts.DispModulo = 1;
    
    % passing everyting into the algorithm
    [xmin,fmin,counteval,stopflag,out] = cmaes('Error',BestGuess,Sigma,opts,Varargin);
end





