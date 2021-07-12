function points = billiard_walk_low_dim(R, r, W, N)
    
    n = length(R);
    L = 2;
    
    A = [ eye(n) ; -eye(n); R];
    b = [ ones(n,1) ; zeros(n,1); r];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];
    
    
    NN = null(Aeq);
    shift = (1/n)*ones(n,1);
    
    b = b - A * shift;
    A = A * NN;
    
    [x, ~]=get_cheb(A,b);
    
    
    sqrt_sum = sqrt(sum(A.^2,2));
    A = A ./ repmat(sqrt_sum, [1 n-1]);
    b = b ./ sqrt_sum;

    d = size(A, 2);
    points = zeros(d, N);

    
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

    points = NN * points + repmat(shift, [1 N]);
end