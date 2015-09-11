% a function which returns the fixed points of the model, given the values
% of the parameters we use . See the mathematica notebook FixedPoint.nb for how these have been
% computed

function x0 = InitialState(theta,Atc,Iptg)

% the list of known parameter values
%degradation parameters
delta_g = 0.0005; %1/min (Andersen et al., Appl. Environ. Microbiol., 1998)
matur = 0.132; %1/min (Iizuka et al., Anal. Chem., 2011)
vz = 100; %nM/min (rate of enzymatic degradation) (Hersch, PNAS, 2004)
Kz = 75; %nM (dissociation constant of enzymatic degradation) (Hersch, PNAS, 2004)
% other parameters
z0 = 9; % experimentally determined

% reading in the parameters we are currently guessing
mu = theta(1);
F = theta(2);
ratio = theta(3);

% set whether atc or iptg have been turned on at a constant level
if Atc==0, fAtc = 0; %aTc
else fAtc = 1;
end;
if Iptg ==0, fIptg = f_lac; %IPTG
else fIptg = 1;
end;

% using the analytical expression for the fixed point, return the state
% values. see mathematica notebook for the expressions
% note this solution makes assumptions on parameter ranges.
%Assumptions -> {am > 0 && dm > 0 && k_on > 0 && k_off > 0 && as > 0 && 
%    ds > 0 && dsm > 0 && dsm > k_off && m > 0 && s > 0 && y > 0}]

% now solve for p and g
% we can get an equation for (p+g) by adding the two eqns
a = F * fAtc;
b = mu + delta_g;
root = roots([b,(vz+Kz*b-a),-Kz*a]);
%just get the positive root
root = root(root>=0);
p = a/(matur + mu + delta_g +(vz/(Kz+root)));
g = root - p;
z = z0 + (g/ratio);
x0 = [p , z ];

% %%a bit of code to check we really are getting the fixed point
% x = x0;
% copies*a_tet/fAtc - mu*x(1) - delta_s*x(1) - k_on*x(1)*x(2) + k_off*x(3) %sRNA
% copies*a_lac/fIptg - mu*x(2) - delta_m*x(2) - k_on*x(1)*x(2) + k_off*x(3) %mRNA
% k_on*x(1)*x(2) - k_off*x(3) - k_hyb*x(3) - mu*x(3) -  delta_m*x(3) %sRNA:mRNA_intermediate
% k_hyb*x(3) - mu*x(4) - (delta_m)*x(4) %sRNA:mRNA_stable
% beta*x(2) + f_srna*beta*x(4) - matur*x(5) - mu*x(5) - delta_g*x(5) - (vz*x(5))/(Kz + x(5) + g) %GFP non-mature
% matur*x(5) - (mu + delta_g)*g - ((vz*g)/(Kz + x(5) + g)) %measured fluoresence
% end





