function XX = hmc_leapfrog(A, b, X0, sigma, mu, a, q, N, step)
    
    d = length(mu);
    Delta_max = 1000;
    
    XX = zeros(d, N);
    
    h = waitbar(0,'Computing samples...');
    for i = 1:N
        
        x_counting_total = 0;
        v = randn(d, 1);
        v_pl = v;
        v_min = -v;
        X_pl = X0;
        X_min = X0;
        %Ti_pl = T0i;
        %Ti_min = T0i;
        %X_nx = X0;
        %Ti_nx = T0i;
        %v_nx = v;
        %grad_x_pl = esti_grad_opt(f_utils, X0, T0i);
        grad_x_pl = a * (2*(sigma*X0) - q*mu); 
        grad_x_min = grad_x_pl;
        %h1 = f(X0, T0i, f_utils) + 0.5 * (v' * v);
        h1 = a * (X0'*sigma*X0 - q*(mu'*X0)) + 0.5 * (v' * v); 
        %exp(-h1)
        uu = log(rand) - h1;
        %x_pl = [X0(lower); T0i];
        %x_min = [X0(lower); T0i];
        %u = rand * exp(-h1);
        j = -1;
        s = 1;
        while(s == 1)
            j = j + 1;
            dir = rand;
            %v = v*dir;
                
            if dir > 0.5
                %v_pl = v_pl - (step/2) * grad_x_pl;
                grad_x = grad_x_pl;
                v = v_pl;
                X = X_pl;
                %Ti = Ti_pl;
            else
                %v_min = v_min - (step/2) * grad_x_min;
                grad_x = grad_x_min;
                v = v_min;
                X = X_min;
                %Ti = Ti_min;
            end
            X_rnd_j = X;
            %Ti_rnd_j = Ti;
            %v_rnd_j = v;
                
            x_counting = 0;
            for k = 1:2^j
                %grad_x
                v = v - (step/2) * grad_x;
                T = step;
                
                lambdas = (b - A*X) ./ A*v;
                [~, pos_max] = min(1./lambdas);
                
                while(true)
                    
                    %v
                    %norm_v = norm(v);
                    %v = v / norm_v;
                    lambdas = (b - A*X) ./ (A*v);
                    %pos_max
                    lambdas(pos_max) = -1;
                    [l_max, pos_max] = max(1./lambdas);
                    l_max = 1 / l_max;
                    
                        
                    % update the current point of the random walk
                    if (T <= l_max)
                        %x = x + T * v;
                        X = X + T * v;
                        %Ti = Ti + T * v((n+1):end);
                        break;
                    end
                        
                    lambda = 0.995*l_max;
                    %x = x + lambda * v;
                    X = X + lambda * v;
                    %Ti = Ti + lambda * v((n+1):end);
                    
                    %reflevt the ray
                    v = v - (2*A(pos_max, :)*v)*A(pos_max, :)';
                    %p = pos(pos_max);
                    %    vv(p) = -vv(p);
                    %    v(1:n) = vv;
                    
                    T = T - lambda;
                end
                
                %grad_x = esti_grad_opt(f_utils, X, Ti);
                grad_x = a * (2*(sigma*X) - q*mu); 
                v = v - (step/2) * grad_x;
                %hj = f(X, Ti, f_utils) + 0.5 * (v' * v);
                hj = a * (X'*sigma*X - q*(mu'*X)) + 0.5 * (v' * v); 
                if (uu > Delta_max - hj)
                    s = 0;
                    break;
                end
                pos_state = false;
                if (uu < -hj)
                    pos_state = true;
                    x_counting = x_counting + 1;
                    x_counting_total = x_counting_total + 1;
                end
                %pos_state
                 
                if (k==1)
                    if (dir > 0.5)
                        X_min_j = X;
                        %Ti_min_j = Ti;
                        v_min_j = v;
                    else
                        X_pl_j = X;
                        %Ti_pl_j = Ti;
                        v_pl_j = v;
                    end
                end
                if (k==2^j)
                    if (dir > 0.5)
                        x_pl = X;
                        x_min = X_min_j;
                        x_pl_min = (x_pl - x_min)';
                        if ((x_pl_min*v < 0) || (x_pl_min*v_min_j < 0))
                            s = 0;
                        end
                    else
                        x_pl = X_pl_j;
                        x_min = X;
                        x_pl_min = (x_pl - x_min)';
                        if ((x_pl_min*v < 0) || (x_pl_min*v_pl_j < 0))
                            s = 0;
                        end
                    end
                end
                if (rand < (1/x_counting) && pos_state)
                    X_rnd_j = X;
                    %Ti_rnd_j = Ti;
                    %v_rnd_j = v;
                end
            end
            if (dir > 0.5)
                X_pl = X;
                v_pl = v;
                %Ti_pl = Ti;
                grad_x_pl = grad_x;
            else
                X_min = X;
                v_min = v;
                %Ti_min = Ti;
                grad_x_min = grad_x;
            end
                
            if (s==1 && rand < (x_counting / x_counting_total))
                X0 = X_rnd_j;
                %T0i = Ti_rnd_j;
                %v_nx = v_rnd_j;
            end
                
            if (s==1)
                x_pl = X_pl;
                x_min = X_min;
                %(x_pl - x_min)*v_min'
                %(x_pl - x_min)*v_pl'
                x_pl_min = (x_pl - x_min)';
                if ((x_pl_min*v_min < 0) || (x_pl_min*v_pl < 0))
                    s = 0;
                end
            end
            waitbar(i/N);
        end
        
        XX(:,i) = X0;
    end
    close(h);
end



