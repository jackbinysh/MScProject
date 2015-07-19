function dx = RibodynamicsModel(t,x, theta)

    % this defines the structure of the vector theta, the list of
    % parameters which goes into the model
    theta = num2cell(theta);
    [copies,a_lac,a_tet,f_srna,f_lac,f_tet,k_on,k_off,k_hyb,delta_m,delta_s,delta_sm,delta_c,delta_g,matur,mu,beta,c,z0,vz,Kz] = theta{:};

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

    dx(1) = copies*a_tet/fu - mu*x(1) - delta_s*x(1) - k_on*x(1)*x(2) + k_off*x(3) + delta_m*(x(3) + x(4)); %sRNA

    dx(2) = copies*a_lac/fv - mu*x(2) - delta_m*x(2) - k_on*x(1)*x(2) + k_off*x(3) + delta_s*(x(3) + x(4)); %mRNA

    dx(3) = k_on*x(1)*x(2) - k_off*x(3) - k_hyb*x(3) - mu*x(3) - (delta_sm + delta_m + delta_s)*x(3); %sRNA:mRNA_intermediate

    dx(4) = k_hyb*x(3) - mu*x(4) - (delta_c + delta_m + delta_s)*x(4); %sRNA:mRNA_stable

    Translation_Rate = beta*x(2) + f_srna*beta*x(4);

    dx(5) = Translation_Rate - matur*x(5) - mu*x(5) - delta_g*x(5) - vz*x(5)/(Kz + x(5) + c*(x(6)-z0)); %GFP non-mature

    dx(6) = (1/c)* (matur*x(5) - mu*(x(6) + delta_g)*c*(x(6)-z0) - (vz*c*(x(6)-z0))/(Kz + x(5) + c*(x(6)-z0))); %measured fluoresence

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


