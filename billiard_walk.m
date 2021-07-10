function [points] = billiard_walk(A, b, x, W, N, L)
    
    d = size(A, 2);
    points = zeros(d, N);
    
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 d]);
    b = b ./ sqrt_sum;    
    
    for i=1:N
        for j=1:W
            
            v = randn(d,1);
            v=v / norm(v);
            T = rand*L;
            while(true)
        
                ls = (b - A*x)./(A*v);
                ls = 1./ls;
                [lmin, facet] = max(ls);
                lmin = 1 / lmin;
                %ls
                %facet = find(ls==lmin);
   
                if (T<=lmin)
                    x = x + T*v;
                    break;
                else
                    x = x + (0.995*lmin)*v;
                    T = T - lmin;
                    v = v -(2*v'*A(facet,:)')*A(facet,:)';
                end
            end
            
    
        end
        points(:,i)=x;
    end

end