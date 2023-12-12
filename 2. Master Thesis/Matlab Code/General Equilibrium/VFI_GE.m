function [v11, v12, v21, v22, ap11, ap12, ap21, ap22, c11, c12, c21, c22] = VFI_GE(r, w1, w2)

% sol the Bellman equation 
% input: equilibrium prices
% output: policy and value functions 

global beta gamma I J P agrid zgrid aa zz gamma11 gamma12 gamma21 gamma22 tol scale_phi

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
c21 = zmat*w1*gamma21 +(1+r)*amat - apmat;  
c22 = zmat*w2*gamma22 +(1+r)*amat - apmat; 

u11 = double(util(c11));
u12 = double(util(c12));
u21 = double(util(c21));
u22 = double(util(c22));

u11(c11 < 0) = - 1e8;
u12(c12 < 0) = - 1e8;
u21(c21 < 0) = - 1e8;
u22(c22 < 0) = - 1e8;

% guess the value function 
v011 = zeros(J,I);
v012 = zeros(J,I);
v021 = zeros(J,I);
v022 = zeros(J,I);

while err > tol && it < 20000

% expected value
Ev = scale_phi.*log(exp(P*v011./scale_phi) + exp(P*v012./scale_phi));
Ev1 = repmat(Ev(1,:),I,1);
Ev2 = repmat(Ev(2,:),I,1);
Ev = [Ev1; Ev2]; 

Ev11 = scale_phi.*log(exp(P*v021./scale_phi) + exp(P*v022./scale_phi));
Ev11 = repmat(Ev1(1,:),I,1);
Ev21 = repmat(Ev1(2,:),I,1);
Ev1 = [Ev11; Ev21]; 

% Bellman operator w\o grid refinement
[v11,ind11] = max(u11 + beta*Ev,[],2);
[v12,ind12] = max(u12 + beta*Ev,[],2);
[v21,ind21] = max(u21 + beta*Ev1,[],2);
[v22,ind22] = max(u22 + beta*Ev1,[],2);

ap11 = agrid(ind11);
ap12 = agrid(ind12);
ap21 = agrid(ind21);
ap22 = agrid(ind22);

% check convergence and update 
v011 = [v011(1,:)'; v011(2,:)'];
v012 = [v012(1,:)'; v012(2,:)'];
v021 = [v021(1,:)'; v021(2,:)'];
v022 = [v022(1,:)'; v022(2,:)'];

err = max(abs(v11 - v011)) + max(abs(v12 - v012)) + max(abs(v21 - v021)) + max(abs(v22 - v022))

v011 = [v11(1:I)'; v11(I+1:I*J)']; 
v012 = [v12(1:I)'; v12(I+1:I*J)']; 
v021 = [v21(1:I)'; v21(I+1:I*J)']; 
v022 = [v22(1:I)'; v22(I+1:I*J)']; 
it = it + 1;  

end

% output
v11 = [v11(1:I) v11(I+1:I*J)]; 
v12 = [v12(1:I) v12(I+1:I*J)]; 
v21 = [v21(1:I) v21(I+1:I*J)]; 
v22 = [v22(1:I) v22(I+1:I*J)]; 

ap11 = [ap11(1:I) ap11(I+1:I*J)]; 
ap12 = [ap12(1:I) ap12(I+1:I*J)]; 
ap21 = [ap21(1:I) ap21(I+1:I*J)]; 
ap22 = [ap22(1:I) ap22(I+1:I*J)]; 

c11 = zz*w1*gamma11 + (1+r)*aa - ap11; 
c12 = zz*w2*gamma12 + (1+r)*aa - ap12;
c21 = zz*w1*gamma21 + (1+r)*aa - ap21; 
c22 = zz*w2*gamma22 + (1+r)*aa - ap22;

fprintf('Value function converged, (it, err) = (%.0f,%.0e) \n',it,err);

end