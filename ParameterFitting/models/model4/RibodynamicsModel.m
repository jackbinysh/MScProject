function dx = RibodynamicsModel(t,x, theta,dataset)

    % the list of known parameter values
    
    %PLtetO1 and PLlacO1 parameters
    f_tet = 2535; %dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    f_lac = 620; %dimensionless (fold) (Lutz and Bujard, NAR, 1997)
    a_tet = 11; %nM/min (Lutz and Bujard, NAR, 1997)
    a_lac = 11; %nM/min (Lutz and Bujard, NAR, 1997)
    
    %degradation parameters

    vz = 100; %nM/min (rate of enzymatic degradation) (Hersch, PNAS, 2004)   
    delta_g = 0.0005; %1/min (Andersen et al., Appl. Environ. Microbiol., 1998)
    matur = 0.132; %1/min (Iizuka et al., Anal. Chem., 2011)
    Kz = 75; %nM (dissociation constant of enzymatic degradation) (Hersch, PNAS, 2004)
    % other parameters
    z0 = 9; % experimentally determined
    copies = 300; %(plasmid copy number)

    % reading in the parameters we are currently guessing
    % order should be {'mu';'F';'c'};
    mu = theta(1);
    F = theta(2);
    c = theta(3);

    % determining the forcing
    t0=0;
    if t<t0
        Atc=0;
        Iptg=0;
    else
        Atc = atc_input(t-t0,dataset); %aTc
        Iptg = iptg_input(t-t0,dataset); %IPTG
    end
    
    % we assume a non zero concentration of atc and iptg saturates the fold
    % so these are just boolean, on/off, either totally repressed or
    % totally unrepressed
    if Atc==0
        fAtc = 0; %aTc
    else
        fAtc = 1;
    end    
    if Iptg == 0
        fIptg = f_lac; %IPTG
    else
        fIptg = 1;
    end;
    
    %%% state equations
 
    dx(1) = F*fAtc - matur*x(1) - mu*x(1) - delta_g*x(1) - (vz*x(1))/(Kz + x(1) + c*(x(2)-z0)); %GFP non-mature

    dx(2) = (1/c)* ( matur*x(1) - (mu + delta_g)*c*(x(2)-z0) - ((vz*c*(x(2)-z0))/(Kz + x(1) + c*(x(2)-z0))) ); %measured fluoresence

    dx = dx';
    
end