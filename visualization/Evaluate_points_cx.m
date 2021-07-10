function [Evals] = Evaluate_points_cx(X, c, a)

N = size(X,2);

Evals = exp(a*(c * X));
   
Evals = Evals/sum(Evals);

end