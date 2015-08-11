function J = Jacobian(t,y,theta,Dataset)

    % the state variables
    s = y(1);
    m = y(2);
    p = y(5);
    z = y(6);
    
    % the list of known parameter values
    %PLtetO1 and PLlacO1 parameters
    f_tet = 2535; %dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    f_lac = 620; %dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    a_tet = 11; %nM/min (Lutz and Bujard, NAR, 1997)
    a_lac = 11; %nM/min (Lutz and Bujard, NAR, 1997)
    %degradation parameters
    delta_g = 0.0005; %1/min (Andersen et al., Appl. Environ. Microbiol., 1998)
    matur = 0.132; %1/min (Iizuka et al., Anal. Chem., 2011)
    vz = 100; %nM/min (rate of enzymatic degradation) (Hersch, PNAS, 2004)
    Kz = 75; %nM (dissociation constant of enzymatic degradation) (Hersch, PNAS, 2004)
    % other parameters
    z0 = 9; % experimentally determined
    copies = 300; %(plasmid copy number)

    % reading in the parameters we are currently guessing
    % order should be 
    %{'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'};
    f_srna = theta(1);
    k_on = theta(2);
    k_off = theta(3);
    k_hyb = theta(4);
    delta_m = theta(5);
    delta_s = theta(6);
    mu = theta(7);
    Beta = theta(8);
    ratio = theta(9);

    % compute the jacobian
    
    J(1,1) = -(mu + delta_s+ k_matur) ;
    J(1,2) = 0;
    J(1,3) = 0;
    J(1,4) = 0;
    J(1,5) = 0;
    J(1,6) = 0;
    J(1,7) = 0;
    
    J(2,1) = k_matur;
    J(2,2) = - (delta_s + mu) - k_on*m
    J(2,3) = -k_on*s
    J(2,4) = k_off;
    J(2,5) = 0;
    J(2,6) = 0;
    J(2,7) =0;
    
    J(3,1) = 0
    J(3,2) = -k_on*m;
    J(3,3) = (mu + delta_m) -k_on*s;
    J(3,4) = k_off;
    J(3,5) = 0;
    J(3,6) = 0;
    J(3,7) = 0;
    
    J(4,1) = 0;
    J(4,2) = k_on*m;
    J(4,3) = k_on*s;
    J(4,4) = -(delta_m+k_off+k_hyb+mu);
    J(4,5) = 0;
    J(4,6) = 0;
    J(4,7) = 0;
    
    J(5,1) = 0;
    J(5,2) = 0;
    J(5,3) = 0
    J(5,4) =  k_hyb;
    J(5,5) = -(mu + delta_m);
    J(5,6) = 0;
    J(5,7) = 0
    
    J(6,1) = 0;
    J(6,2) = 0;
    J(6,3) = Beta;
    J(6,4) = 0;
    J(6,5) = f_srna * Beta;
    J(6,6) = -(matur + mu +delta_g) -(vz/(Kz + ratio*(z-z0) + p)) + (vz*p)*(1/(Kz + ratio*(z-z0) + p))^2;
    J(6,7) = ratio*( (vz*p)*(1/(Kz + ratio*(z-z0) + p))^2 );
    
    J(7,1) = 0;
    J(7,2) = 0;
    J(7,3) = 0;
    J(7,4) = 0;
    J(7,5) = 0;
    J(7,6) = (1/ratio) * ( matur  + (vz*ratio*(z-z0))*(1/(Kz + ratio*(z-z0) + p))^2 );
    J(7,7) = (1/ratio) *  ( - ratio*(mu + delta_g) - (vz*ratio/(Kz + ratio*(z-z0) + p)) + (vz*(ratio^2)*(z-z0))*(1/(Kz + ratio*(z-z0) + p))^2 );
    

end

