function [vals2] = Evaluate_mixture3(X, S, mu, a, q, w)

M = length(w);
N = size(X,2);

Evals =  mu * X;
vals = zeros(M,N);

Sx = S*X;



   
    for j=1:M
        for i=1:N
           xS = X(:,i)' * Sx(:,i);
           vals(j,i) = exp( -a(j) * ( -q(j)*Evals(j,i) + xS ) );
        end
        vals(j,:) = vals(j,:) / (mean(vals(j,:))*0.5);
    end
    
    vals2 = zeros(1,N);
    for i=1:N
        for j=1:M
          vals2(i) = vals2(i) + w(j) * vals(j,i);
        end
    end


end