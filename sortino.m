function [S,DD] = sortino(R, x0, sigma, mu)
%   [S,DD] = sortino(R, MAR) returns the Sortino Ratio and Downside
%   Deviation using the historical returns in vector R and a minimum
%   acceptable return (MAR)
%  ======================================================================
%  INPUTS:
%   R    - vector of fund returns
%   MAR  - minimum acceptable return
%   
%  OUTPUTS:
%   DD   - downside deviation
%   S    - sortino ratio
%  ======================================================================
%
%   Author: Lorenzo Brancali
%   E-mail: lbrancali@gmail.com
%   Date:   5th March 2012
%
%  ======================================================================

n = length(mu);
EqPtf = ones(n,1)/n;

A = [ eye(n) ; -eye(n)];
b = [ ones(n,1) ; zeros(n,1)];
Aeq = [ ones(1,n) ];
beq = [ 1 ];
    
options = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);
    
Ptf = fmincon(@(x) volatility_fun(x, sigma), EqPtf, A, b, Aeq, beq, [], [], [], options);
MAR = mu'*Ptf;
ret = R'*x0;
for j=1:size(R,2) 
    
    F=R(:,j);
    DD(:,j)=sqrt(sum(nonzeros(F(F<MAR)).^2))/length(nonzeros(F(F<MAR)));
    S(:,j)=(ret-MAR)/DD(:,j);
end