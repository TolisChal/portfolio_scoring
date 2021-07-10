function [Evals] = Evaluate_points(X, S, mu, sigma)

N = size(X,2);

Evals = mu * X;

Sx = S*X;

for i=1:N
   
    Evals(i) = exp(sigma * (Evals(i) / (X(:,i)' * Sx(:,i))));
    
end
   
Evals = Evals/sum(Evals);

end