function faz = CK_GE(ap)

% iterate over Chapman-Kolmogorov equation to find the stationary 
% distribution for (a,z)

global I J P agrid tol

% construct the transition matrix from each (a,z) to each (a',z')
T = zeros(I*J,I*J); 
Ia = zeros(I*J,I);
PP = kron(P,ones(I,1));
% To find a' on the asset grid. Given j split the probability 1 of the move from a to a' to
% the two agrid values (a-1, a+1) closest to ap so that (1-p)a-1 + pa1 = a'. This is useful 
% for EGM, in our case by construction we are always on the grid when moving from a to a' 
policy = ap(:); 
for s = 1:I*J
% find a' on the asset grid   
[val, ind] = position(agrid,I,policy(s));   
Ia(s,:) = val(1)*(agrid == agrid(ind(1))) + val(2)*(agrid == agrid(ind(2)));    
% find the transitions probabilities
T(s,:) = kron(PP(s,:),Ia(s,:));  
end
T = sparse(T); 

% find stationary distribution by fixed point iteration 
faz = (1/(I*J))*ones(I*J,1);
err = 1; 
while err > tol
    faz1 = T'*faz;
    err = max(abs(faz1-faz));
    faz = faz1;
end

faz = reshape(faz,I,J); 

fprintf('Stationary distribution found \n');

end


function [val, ind] = position(grid,n,x)

% returns the values and the indices of the two points on column vector grid, 
% with dimension n and ascending order, that are closets to the scalar x 
% if x is on the grid return only one nonzero value

if min(abs(grid - x)) == 0 % x is on the grid already

ind(1) = find(grid == x);
ind(2) = find(grid == x);
val(1) = 1;
val(2) = 0;    
    
else % x is not on the grid 

% find lower bound on the grid
[~,ind] = sort([grid; x]);
temp = find(ind > n);
i = (temp - 1);

if ((i+1)>n)
ind(1) = n;
ind(2) = n;
val(1) = 1;
val(2) = 0;
elseif (i==0)
ind(1) = 1;
ind(2) = 1;
val(1) = 1;
val(2) = 0;
else
ind(1) = i;
ind(2) = i+1;
dist = grid(i+1)-grid(i);
val(2)=(x - grid(i))/dist;
val(1)=(grid(i+1)- x)/dist;
end

end
end


