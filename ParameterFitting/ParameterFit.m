function [xmin,fmin,counteval,stopflag,out] = ParameterFit(theta,x0,Dataset)

    % reading off the initial guess and the bounds
    BestGuess = [x0. BestGuess; theta.BestGuess];
    opts.LBounds = [x0.LowerBounds; theta.LowerBounds];  
    opts.UBounds = [x0.UpperBounds; theta.UpperBounds];
    
    % The cmaes wants all parameters to be on the same scale. I take our
    % initial Guess as setting the scale for the problem
    Scale =  BestGuess;
    opts.LBounds = opts.LBounds./Scale;
    opts.UBounds = opts.UBounds./Scale;
    BestGuess = BestGuess./Scale;
    
    % setting an initial spread of points, using 1/3rd of the search region as recommended in cmaes script
    Sigma = (opts.UBounds - opts.LBounds)/3; 
    
    % setting some display etc. options
    opts.DispModulo = 1;
    % constructing the struct of additional arguments to pass to the
    % objective
    Varargin = struct('Scale',Scale,'Dataset',Dataset);
    [xmin,fmin,counteval,stopflag,out] = cmaes('Error',BestGuess,Sigma,opts,Varargin);
end





