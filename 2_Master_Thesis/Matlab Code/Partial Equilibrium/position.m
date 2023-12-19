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
