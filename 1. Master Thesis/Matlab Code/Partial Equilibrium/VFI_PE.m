function [v11, v12, ap11, ap12, c11, c12] = VFI_PE(r,w1, w2)

% sol the Bellman equation 
% input: equilibrium prices
% output: policy and value functions 

global beta gamma I J P agrid zgrid aa zz gamma11 gamma12 scale_phi

% initialization 
it = 0;
err = 1; 
apmat = repmat(agrid',I*J,1); 
amat = repmat(agrid,J,I); 
zmat = kron(zgrid,ones(I,I));

% utility functions 
if gamma > 1
util = @(x) (x.^(1-gamma))./(1-gamma); 
elseif gamma == 1
util = @(x) log(x); 
end

% compute utility IJ x I matrix 
c11 = zmat*w1*gamma11 +(1+r)*amat - apmat;  
c12 = zmat*w2*gamma12 +(1+r)*amat - apmat; 

u11 = double(util(c11));
u12 = double(util(c12));

u11(c11 < 0) = - 1e8;
u12(c12 < 0) = - 1e8;

% guess the value function 
v011 = zeros(J,I);
v012 = zeros(J,I);

while err > 1e-13 && it < 20000

% expected value
Ev = scale_phi.*log(exp(P*v011./scale_phi) + exp(P*v012./scale_phi));
Ev1 = repmat(Ev(1,:),I,1);
Ev2 = repmat(Ev(2,:),I,1);
Ev = [Ev1; Ev2]; 

% Bellman operator w\o grid refinement
[v11,ind11] = max(u11 + beta*Ev,[],2);
[v12,ind12] = max(u12 + beta*Ev,[],2);

ap11 = agrid(ind11);
ap12 = agrid(ind12);

% check convergence and update 
v011 = [v011(1,:)'; v011(2,:)'];
v012 = [v012(1,:)'; v012(2,:)'];
err = max(abs(v11 - v011)) + max(abs(v12 - v012))
v011 = [v11(1:I)'; v11(I+1:I*J)']; 
v012 = [v12(1:I)'; v12(I+1:I*J)']; 
it = it + 1;  

end

% output
v11 = [v11(1:I) v11(I+1:I*J)]; 
v12 = [v12(1:I) v12(I+1:I*J)]; 

ap11 = [ap11(1:I) ap11(I+1:I*J)]; 
ap12 = [ap12(1:I) ap12(I+1:I*J)]; 

c11 = zz*w1*gamma11 + (1+r)*aa - ap11; 
c12 = zz*w2*gamma12 + (1+r)*aa - ap12;

fprintf('Value function converged, (it, err) = (%.0f,%.0e) \n',it,err);

end