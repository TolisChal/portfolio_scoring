function [vals] = Evaluate_mixture1(X, S, mu, a, w)

N = size(X,2);
M = length(w);

Evals = mu * X;
vals = zeros(1,N);

Sx = S*X;

for i=1:N
    
    xS = X(:,i)' * Sx(:,i);
    for j=1:M
        vals(i) = vals(i) + w(j) * exp( a(j) * ( Evals(i) / xS ) );
    end
    
end
   
vals = vals/sum(vals);

end