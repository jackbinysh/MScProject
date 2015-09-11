function [List,pset] = Identifiability (SSM)
%this function uses SSM to to determine identifiability of parameters

[r,c]=size(SSM);
X=[];
pset=[];

P_i_name = {'f_srna';'k_on';'k_off';'k_hyb';'delta_m';'delta_s';'mu';'beta';'c'};

% McAuley procedure doi:10.1081/PRE-120024426
R=SSM; % the first time, use SSM to do column sumComputeIdentifiability

for j=1:c
     for i=1:c
         M(i)=R(:,i)'*R(:,i); % the square sum of each column
     end
     
     [a,pos]=max(M);          % finds the indices of the maximum values of M, and returns them in output vector pos. 
     
     if a>1.0e-07                 % this is the tolerance
         X=[X SSM(:,pos)];      % the colomn that has the largest SS magnitude
         pset=[pset; pos];      % give the index of the parameter
         Shat=X*inv(X'*X)*X'*SSM; % Find the prediction SSM
         R=SSM-Shat;            % residual mtrx, now the residual matrix is the new mtrx that we find the next identifiable parameter, return to the top of the J loop  
     end
 end

pset=unique(pset,'stable');                % identifiable parameters

c=size(pset);
for i=1:c
    k=pset(i);
    List(i,:)=P_i_name(k,:);
end

end

