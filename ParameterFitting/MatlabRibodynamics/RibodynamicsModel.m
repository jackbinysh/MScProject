function dx = RibodynamicsModel(t,x, theta)

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
delta_sm = 0; %1/min

% reading in the parameters we are currently guessing
f_srna = theta{'f_srna',1};
k_on = theta{'k_on',1};
k_off = theta{'k_off',1};
k_hyb = theta{'k_hyb',1};
delta_m = theta{'delta_m',1};
delta_s = theta{'delta_s',1};
delta_c = theta{'delta_c',1};
mu = theta{'mu',1};
beta = theta{'beta',1};
c = theta{'c',1};

    %%% determining the forcing
    t0=0;
    if t<t0,
        u=0;
        v=0;
    else
        u = atc_input(t-t0); %aTc
        v = iptg_input(t-t0); %IPTG
    end;
    if u==0, fu = f_tet; %aTc
    else fu = 1;
    end;
    if v ==0, fv = f_lac; %IPTG
    else fv = 1;
    end;
    %%% state equations

    dx(1) = copies*a_tet/fu - mu*x(1) - delta_s*x(1) - k_on*x(1)*x(2) + k_off*x(3); %sRNA

    dx(2) = copies*a_lac/fv - mu*x(2) - delta_m*x(2) - k_on*x(1)*x(2) + k_off*x(3); %mRNA

    dx(3) = k_on*x(1)*x(2) - k_off*x(3) - k_hyb*x(3) - mu*x(3) - (delta_sm)*x(3); %sRNA:mRNA_intermediate

    dx(4) = k_hyb*x(3) - mu*x(4) - (delta_c)*x(4); %sRNA:mRNA_stable

    Translation_Rate = beta*x(2) + f_srna*beta*x(4);

    dx(5) = Translation_Rate - matur*x(5) - mu*x(5) - delta_g*x(5) - vz*x(5)/(Kz + x(5) + c*(x(6)-z0)); %GFP non-mature

    dx(6) = (1/c)* ( matur*x(5) - (mu + delta_g)*c*(x(6)-z0) - ((vz*c*(x(6)-z0))/(Kz + x(5) + c*(x(6)-z0))) ); %measured fluoresence

    dx = dx';
    
end

    function u = atc_input(t)

        p1=120; %min 
        p2=120; %min
        p=p1+p2;
        tp = t-floor(t/p)*p;
        if tp <= p1,
            u = 2; %ng/mL
        else
            u = 0; %ng/mL
        end
    end

    function v = iptg_input(t)
        p1=120; %min 
        p2=120; %min
        p=p1+p2;
        tp = t-floor(t/p)*p;
        if tp <= p1,
            v = 1; %mM
        else
            v = 1; %mM
        end
    end


