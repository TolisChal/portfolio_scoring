function [vals] = Evaluate_mixture2(X, S, mu, a, q, w)

M = length(w);
N = size(X,2);

Evals =  mu * X;
vals = zeros(1,N);

Sx = S*X;

for i=1:N

    xS = X(:,i)' * Sx(:,i);
    for j=1:M
        vals(i) = vals(i) + w(j) * exp( -a(j) * ( -q(j)*Evals(i) + xS ) );
    end
    
end

vals = vals/sum(vals);

end