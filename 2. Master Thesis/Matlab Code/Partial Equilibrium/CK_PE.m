function faz = CK_PE(ap)

% iterate over Chapman-Kolmogorov equation to find the stationary 
% distribution for (a,z)

global I J P agrid 

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
    [val, ind] = position(agrid, I, policy(s));   
    Ia(s,:) = val(1)*(agrid == agrid(ind(1))) + val(2)*(agrid == agrid(ind(2)));    
    % find the transitions probabilities
    T(s,:) = kron(PP(s,:),Ia(s,:));  
end
T = sparse(T); 

% find stationary distribution by fixed point iteration 
faz = (1/(I*J))*ones(I*J,1);
err = 1; 
while err > 1e-16
    faz1 = T'*faz;
    err = max(abs(faz1-faz))
    faz = faz1;
end

faz = reshape(faz,I,J); 

fprintf('Stationary distribution found \n');

end


