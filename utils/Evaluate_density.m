function [vals] = Evaluate_density(X, sigma, mu, a, q)

N = size(X,2);

Evals =  mu' * X;
vals = zeros(1,N);

Sx = sigma*X;

for i=1:N
    
    xS = X(:,i)' * Sx(:,i);
    vals(i) = exp( -a * ( -q*Evals(i) + xS ) );
end

end