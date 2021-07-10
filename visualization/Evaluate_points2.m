function [Evals] = Evaluate_points2(X, S, mu, sigma)

N = size(X,2);

Evals = mu * X;

Sx = S*X;

for i=1:N
   
    Evals(i) = sigma * (Evals(i) - (X(:,i)' * Sx(:,i)));
    
end
   


end