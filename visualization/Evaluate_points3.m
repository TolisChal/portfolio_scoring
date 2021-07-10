function [Evals] = Evaluate_points3(X, S, mu, sigma, q)

N = size(X,2);

Evals = - q * (mu * X);

Sx = S*X;

for i=1:N
   
    Evals(i) = exp(-sigma * (Evals(i) + (X(:,i)' * Sx(:,i))));
    
end
   
Evals = Evals/sum(Evals);

end