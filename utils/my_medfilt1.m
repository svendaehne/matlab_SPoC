function [x_flt, idx_flt] = my_medfilt1(x, n)

if nargin < 2
    n = 3;
end

x = x(:);

nx = length(x);
if rem(n,2)~=1    % n even
    m = n/2;
else
    m = (n-1)/2;
end

X = [zeros(m,1); x; zeros(m,1)];

indr = (0:n-1)';
indc = 1:nx;

ind = indc(ones(1,n),1:nx) + ...
    indr(:,ones(1,nx));
X = X(ind);
x_flt = median(X)';

[foo, idx] = min(abs(X - repmat(x_flt',[n, 1])));
idx_flt = ind(idx + (0:n:(nx-1)*n)) - m;









