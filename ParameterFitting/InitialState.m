% a function which returns the fixed points of the model, given the values
% of the parameters we use
% see the mathematica notebook FixedPoint.nb for how these have been
% computed

function x0 = InitialState(theta)

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
f_srna = theta(1);
k_on = theta(2);
k_off = theta(3);
k_hyb = theta(4);
delta_m = theta(5);
delta_s = theta(6);
delta_c = theta(7);
mu = theta(8);
beta = theta(9);
ratio = theta(10);

% using the analytical expression for the fixed point, return the state
% values. see mathematica notebook for the expressions
% note this solution makes assumptions on parameter ranges.
%Assumptions -> {am > 0 && dm > 0 && k_on > 0 && k_off > 0 && as > 0 && 
%    ds > 0 && dsm > 0 && dsm > k_off && m > 0 && s > 0 && y > 0}]

am = copies*a_lac/f_lac;
dm = mu + delta_m;
as = copies*a_tet/f_tet;
ds = mu + delta_s;
dsm = k_off + k_hyb + mu +delta_m + delta_s;

m = (1/(2*dm*(dsm-k_off)*k_on))*( -dm*ds*dsm +(am-as)*(dsm-k_off)*k_on + sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-k_off)*k_on + (am-as)^2*(dsm -k_off)^2*k_on^2) );
s = (1/(2*ds*(dsm-k_off)*k_on))*( -dm*ds*dsm - (am-as)*(dsm-k_off)*k_on + sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-k_off)*k_on + (am-as)^2*(dsm -k_off)^2*k_on^2) );
y = (1/(2*(dsm-k_off)^2*k_on))*( dm*ds*dsm +(am+as)*(dsm-k_off)*k_on - sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-k_off)*k_on + (am-as)^2*(dsm -k_off)^2*k_on^2) );
c = (k_hyb/(mu + delta_c))*y;

% now solve for p and g
% we can get an equation for (p+g) by adding the two eqns
a = beta*m + f_srna*beta*c;
b = mu + delta_g;
root = roots([b,(vz+Kz*b-a),-Kz*a]);
%just get the positive root
root = root(root>0);
p = a/(matur + mu + delta_g +(vz/(Kz+root)));
g = root - p;
z = z0 + (g/ratio);
x0 = [s, m ,y, c , p , z ];

% a bit of code to check we really are getting the fixed point
%x = x0;
%copies*a_tet/f_tet - mu*x(1) - delta_s*x(1) - k_on*x(1)*x(2) + k_off*x(3) %sRNA
%copies*a_lac/f_lac - mu*x(2) - delta_m*x(2) - k_on*x(1)*x(2) + k_off*x(3) %mRNA
%k_on*x(1)*x(2) - k_off*x(3) - k_hyb*x(3) - mu*x(3) - (delta_s + delta_m)*x(3) %sRNA:mRNA_intermediate
%k_hyb*x(3) - mu*x(4) - (delta_c)*x(4) %sRNA:mRNA_stable
%beta*x(2) + f_srna*beta*x(4) - matur*x(5) - mu*x(5) - delta_g*x(5) - (vz*x(5))/(Kz + x(5) + g) %GFP non-mature
%matur*x(5) - (mu + delta_g)*g - ((vz*g)/(Kz + x(5) + g)) %measured fluoresence
end





