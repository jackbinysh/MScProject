% a function which returns the fixed points of the model, given the values
% of the parameters we use

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
delta_sm = 0; %1/min

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
%Assumptions -> {am > 0 && dm > 0 && kon > 0 && koff > 0 && as > 0 && 
%    ds > 0 && dsm > 0 && dsm > koff && m > 0 && s > 0 && y > 0}]

am = N*a_lac/f_lac;
dm = mu + delta_m
as = N*a_tet/f_tet;
ds = mu + delta_s;
dsm = k_off + k_on + mu +delta_m + delta_s;

m = (-1/(2*dm*(dsm-koff)*kon))*( (dm*ds*dsm -(am-as)*(dsm-koff)*kon + sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-koff)*kon + (am-as)^2*(dsm -koff)^2*kon^2]) );
s = (-1/(2*ds*(dsm-koff)*kon))*( (dm*ds*dsm +(am-as)*(dsm-koff)*kon + sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-koff)*kon + (am-as)^2*(dsm -koff)^2*kon^2]) );
y = (-1/(2*(dsm-koff)^2*kon))*( (dm*ds*dsm +(am+as)*(dsm-koff)*kon + sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-koff)*kon + (am-as)^2*(dsm -koff)^2*kon^2]) );
c = (k_hyb/(mu + delta_c))*y0;

% now solve for p and g
% we can get an equation for (p+g) by adding the two eqns
a = beta*m - f_s*beta*c;
b = mu + delta_g;
x = roots([b,(vz+Kz*b-a),-Kz*a]);
p = a/(matur + mu + delta_g +(vz/Kz+x));
g = x - p;
z = z0 + (g/ratio);

x0 = [s, m ,y, c , p , z ];




