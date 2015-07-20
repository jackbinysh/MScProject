function [History] = ParameterFit(InitialGuess,UpperBound,LowerBound)

    History.x = [];
    History.fval = [];
    options = optimset('Display','iter','OutputFcn',@outfun, 'MaxIter', 2000);
    [x,fval,exitflag,output] = fminsearch(@Error,InitialGuess, options);
    %[x,fval,exitflag,output] = fmincon(@Error,InitialGuess,[],[],[],[],lb,ub)
    %[x,fval,exitflag,output] = ga(@Error,length(lb),[],[],[],[],lb,ub);

    %% output function for the optimiser
    function stop = outfun(x,optimValues,state)
         stop = false;
         if state == 'iter'
             % Concatenate current point and objective function
             % value with history. x must be a row vector.
               History.fval = [History.fval; optimValues.fval];
               History.x = [History.x; x];
         end
    end

end





