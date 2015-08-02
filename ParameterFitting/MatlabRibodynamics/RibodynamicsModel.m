function dx = RibodynamicsModel(t,x, theta,dataset)

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
    %{'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'delta_c';'mu';'beta';'c'};
    f_srna = theta(1);
    k_on = theta(2);
    k_off = theta(3);
    k_hyb = theta(4);
    delta_m = theta(5);
    delta_s = theta(6);
    delta_c = theta(7);
    mu = theta(8);
    beta = theta(9);
    c = theta(10);

    %%% determining the forcing
    t0=0;
    if t<t0,
        u=0;
        v=0;
    else
        u = atc_input(t-t0,dataset); %aTc
        v = iptg_input(t-t0,dataset); %IPTG
    end
    
    % we assume a non zero concentration of atc and iptg saturates the fold
    % so these are just boolean, on/off, either totally repressed or
    % totally unrepressed
    if u==0, fu = f_tet; %aTc
    else fu = 1;
    end;
    if v ==0, fv = f_lac; %IPTG
    else fv = 1;
    end;
    
    %%% state equations

    dx(1) = copies*a_tet/fu - mu*x(1) - delta_s*x(1) - k_on*x(1)*x(2) + k_off*x(3); %sRNA

    dx(2) = copies*a_lac/fv - mu*x(2) - delta_m*x(2) - k_on*x(1)*x(2) + k_off*x(3); %mRNA

    dx(3) = k_on*x(1)*x(2) - k_off*x(3) - k_hyb*x(3) - mu*x(3) - (delta_s + delta_m)*x(3); %sRNA:mRNA_intermediate

    dx(4) = k_hyb*x(3) - mu*x(4) - (delta_c)*x(4); %sRNA:mRNA_stable
 
    dx(5) = beta*x(2) + f_srna*beta*x(4) - matur*x(5) - mu*x(5) - delta_g*x(5) - (vz*x(5))/(Kz + x(5) + c*(x(6)-z0)); %GFP non-mature

    dx(6) = (1/c)* ( matur*x(5) - (mu + delta_g)*c*(x(6)-z0) - ((vz*c*(x(6)-z0))/(Kz + x(5) + c*(x(6)-z0))) ); %measured fluoresence

    dx = dx';
    
end