% a function which returns the fixed points of the model, given the values
% of the parameters we use . See the mathematica notebook FixedPoint.nb for how these have been
% computed

function x0 = InitialState(theta,Atc,Iptg)

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
 % order should be {'f_srna';'k_on';'k_off';'k_hyb';'k_matur';'delta_m';'delta_s';'mu';'beta';'c'};
f_srna = theta(1);
k_on = theta(2);
k_off = theta(3);
k_hyb = theta(4);
k_matur = theta(5);
delta_m = theta(6);
delta_s = theta(7);
mu = theta(8);
beta = theta(9);
ratio = theta(10);


% set whether atc or iptg have been turned on at a constant level
if Atc==0, fAtc = f_tet; %aTc
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

si = ( copies*a_tet/fAtc ) / ( mu+delta_s + k_matur );

am = copies*a_lac/fIptg;
dm = mu + delta_m;
as = k_matur*si;
ds = mu + delta_s;
dsm = k_off + k_hyb + mu +delta_m;

m = (1/(2*dm*(dsm-k_off)*k_on))*( -dm*ds*dsm +(am-as)*(dsm-k_off)*k_on + sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-k_off)*k_on + (am-as)^2*(dsm -k_off)^2*k_on^2) );
s = (1/(2*ds*(dsm-k_off)*k_on))*( -dm*ds*dsm - (am-as)*(dsm-k_off)*k_on + sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-k_off)*k_on + (am-as)^2*(dsm -k_off)^2*k_on^2) );
y = (1/(2*(dsm-k_off)^2*k_on))*( dm*ds*dsm +(am+as)*(dsm-k_off)*k_on - sqrt(dm^2*ds^2*dsm^2 + 2*(am+as)*dm*ds*dsm*(dsm-k_off)*k_on + (am-as)^2*(dsm -k_off)^2*k_on^2) );
c = (k_hyb/(mu + delta_m))*y;

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
x0 = [si, s, m ,y, c , p , z ];


% %a bit of code to check we are getting the fixed point 
% x = x0;
% 
%     dx(1) = copies*a_tet/fAtc - (mu + delta_s)*x(1) - k_matur*x(1); %immature sRNA
% 
%     dx(2) = k_matur*x(1) - (mu+ delta_s) *x(2) - k_on*x(2)*x(3) + k_off*x(4); %mature sRNA
% 
%     dx(3) = copies*a_lac/fIptg - (mu + delta_m)*x(3)- k_on*x(2)*x(3) + k_off*x(4); %mRNA
%     
%     dx(4) = k_on*x(2)*x(3) - k_off*x(4) - k_hyb*x(4) - mu*x(4) - delta_m*x(4); %sRNA:mRNA_intermediate
%     
%     dx(5) = k_hyb*x(4) - mu*x(5) - (delta_m)*x(5); %sRNA:mRNA_stable
%      
%     dx(6) = beta*x(3) + f_srna*beta*x(5) - matur*x(6) - mu*x(6) - delta_g*x(6) - (vz*x(6))/(Kz + x(6) + ratio*(x(7)-z0)); %GFP non-mature
% 
%     dx(7) = (1/ratio)* ( matur*x(6) - (mu + delta_g)*ratio*(x(7)-z0) - ((vz*ratio*(x(7)-z0))/(Kz + x(6) + ratio*(x(7)-z0))) ); %measured fluoresence
%  
%     dx = dx';

end





