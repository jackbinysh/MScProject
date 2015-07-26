function ribodynamics()

dim=6;
tf = 700; %min
y0 = zeros(dim,1);
[t,y] = ode23s(@model,[0 tf],y0,[]);
z0=9; %autofluorescence;
c=973; %ratio concentration/fluorescence TO BE FITTED
z=z0+y(:,6)/c;
atc = zeros(length(t),1);
iptg = zeros(length(t),1);
for i=1:length(t),
    atc(i)=atc_input(t(i));
    iptg(i)=iptg_input(t(i));
end;
plot(t,z,'k','LineWidth',3);

hold on;
s_tot = y(:,1)+y(:,3)+y(:,4);
m_tot = y(:,2)+y(:,3)+y(:,4);
plot(t,s_tot,'b','LineWidth',3);
plot(t,m_tot,'r','LineWidth',3);
plot(t,y(:,4)./m_tot*10,'r','LineWidth',3);

plot(t,atc+8,'r--','LineWidth',3);
plot(t,iptg,'k--','LineWidth',3);

function dy = model(t,y)

%riboregulation parameters
k_on = 1e5; %1/(nM*min) TO BE FITTED
k_off = 1e7; %1/min TO BE FITTED 
k_hyb = 0.01; %1/min TO BE FITTED
f_srna = 435; %dimensionless (fold) TO BE FITTED

%degradation parameters
mu = 0.0166; %1/min TO BE FITTED
delta_s = 0.1; %1/min TO BE FITTED
delta_m = 0.175; %1/min TO BE FITTED
delta_c = 0.1; %1/min TO BE FITTED
delta_sm = 0; %1/min
delta_g = 0.0005; %1/min (Andersen et al., Appl. Environ. Microbiol., 1998)
matur = 0.132; %1/min (Iizuka et al., Anal. Chem., 2011)
vz = 100; %nM/min (rate of enzymatic degradation) (Hersch, PNAS, 2004)
Kz = 75; %nM (dissociation constant of enzymatic degradation) (Hersch, PNAS, 2004)

%PLtetO1 and PLlacO1 parameters
f_tet = 2535; %dimensionless (fold) (Lutz and Bujard, NAR, 1997)
f_lac = 620; %dimensionless (fold) (Lutz and Bujard, NAR, 1997)
a_tet = 11; %nM/min (Lutz and Bujard, NAR, 1997)
a_lac = 11; %nM/min (Lutz and Bujard, NAR, 1997)

%other parameters
copies = 300; %(plasmid copy number)
beta = 0.0023; %nM/min (translation rate) TO BE FITTED

t0=0;
if t<t0,
    u=0;
    v=0;
else
    u = atc_input(t-t0); %aTc
    v = iptg_input(t-t0); %IPTG
end;

if u==0, fu = f_tet;
else fu = 1;
end;
if v==0, fv = f_lac;
else fv = 1;
end;

dy = zeros(6,1);

dy(1) = copies*a_tet/fu - mu*y(1) - delta_s*y(1) - k_on*y(1)*y(2) + k_off*y(3) + delta_m*(y(3) + y(4)); %sRNA

dy(2) = copies*a_lac/fv - mu*y(2) - delta_m*y(2) - k_on*y(1)*y(2) + k_off*y(3) + delta_s*(y(3) + y(4)); %mRNA

dy(3) = k_on*y(1)*y(2) - k_off*y(3) - k_hyb*y(3) - mu*y(3) - (delta_sm + delta_m + delta_s)*y(3); %sRNA:mRNA_intermediate

dy(4) = k_hyb*y(3) - mu*y(4) - (delta_c + delta_m + delta_s)*y(4); %sRNA:mRNA_stable

Translation_Rate = beta*y(2) + f_srna*beta*y(4);

dy(5) = Translation_Rate - matur*y(5) - mu*y(5) - delta_g*y(5) - vz*y(5)/(Kz + y(5) + y(6)); %GFP non-mature

dy(6) = matur*y(5) - mu*y(6) - delta_g*y(6) - vz*y(6)/(Kz + y(5) + y(6)); %GFP mature

function u = atc_input(t)

p1=120; %min 
p2=120; %min
p=p1+p2;
tp = t-floor(t/p)*p;
if tp <= p1,
    u = 2; %ng/mL
else
    u = 0; %ng/mL
end;

function v = iptg_input(t)

p1=120; %min 
p2=120; %min
p=p1+p2;
tp = t-floor(t/p)*p;
if tp <= p1,
    v = 1; %mM
else
    v = 1; %mM
end;
