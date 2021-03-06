function [c,r]=get_cheb(A,b)

    
    %CHEBYCENTERR Compute Chebyshev center of polytope Ax <= b.
%  r0 acceptable threshold (lower) for r. Optional.
%  The Chebyshev center of a polytope is the center of the largest
%  hypersphere enclosed by the polytope. 
%  Requires optimization toolbox.
[n,p] = size(A);
an = sqrt(sum(A.^2,2));
A1 = zeros(n,p+1);
A1(:,1:p) = A;
A1(:,p+1) = an;
f = zeros(p+1,1);
f(p+1) = -1;
%options = optimset;
%options = optimset(options,'Display', 'off');
%if nargin == 3
    lb = ones(p+1,1) * -Inf;
    ub = ones(p+1,1) * Inf;
    lb(p+1) = 0; ub(p+1) = inf;
%else
%    lb = []; ub = [];
%end
%options = mskoptimset('Display', 'off');
%c = lin_mosek2(f,A1,b,[],[],lb,ub,[]);
options = optimoptions('linprog','Display','off');
c = linprog(f,A1,b,[],[],lb,ub,options);
r = c(p+1);
c = c(1:p);
end