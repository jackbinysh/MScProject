function J = Jacobian(t,y,theta,Dataset)

    % the state variables
    p = y(1);
    z = y(2);
    
    % the list of known parameter values
    %degradation parameters
    delta_g = 0.0005; %1/min (Andersen et al., Appl. Environ. Microbiol., 1998)
    matur = 0.132; %1/min (Iizuka et al., Anal. Chem., 2011)
    vz = 100; %nM/min (rate of enzymatic degradation) (Hersch, PNAS, 2004)
    Kz = 75; %nM (dissociation constant of enzymatic degradation) (Hersch, PNAS, 2004)
    % other parameters
    z0 = 9; % experimentally determined

    % reading in the parameters we are currently guessing
    % order should be 
    %{'mu';'F';'c'};

    mu = theta(1);
    ratio = theta(3);

    % compute the jacobian
        
    J(1,1) = -(matur + mu +delta_g) -(vz/(Kz + ratio*(z-z0) + p)) + (vz*p)*(1/(Kz + ratio*(z-z0) + p))^2;
    J(1,2) = ratio*( (vz*p)*(1/(Kz + ratio*(z-z0) + p))^2 );
    
    J(2,1) = (1/ratio) * ( matur  + (vz*ratio*(z-z0))*(1/(Kz + ratio*(z-z0) + p))^2 );
    J(2,2) = (1/ratio) *  ( - ratio*(mu + delta_g) - (vz*ratio/(Kz + ratio*(z-z0) + p)) + (vz*(ratio^2)*(z-z0))*(1/(Kz + ratio*(z-z0) + p))^2 );
    

end

