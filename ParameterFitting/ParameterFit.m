function [History,output] = ParameterFit(theta,x0)

    History.x = [];
    History.fval = [];

    LowerBounds = [x0.LowerBounds; theta.LowerBounds]';
    BestGuess = [x0. BestGuess; theta.BestGuess]';
    UpperBounds = [x0.UpperBounds; theta.UpperBounds]';

    %options = gaoptimset('Display','iter','OutputFcns',@outfun, 'InitialPopulation',BestGuess);
    options = gaoptimset('Display','iter', 'InitialPopulation',BestGuess);
    
    [x,fval,exitflag,output] = ga(@Error,length(BestGuess),[],[],[],[],...
    LowerBounds,UpperBounds,[],options);

    %% output function for the optimiser
    function state = outfun(x,optimValues,state)
         stop = false;
         if state == 'iter'
             % Concatenate current point and objective function
             % value with history. x must be a row vector.
               History.fval = [History.fval; optimValues.fval];
               History.x = [History.x; x];
         end
    end

end





