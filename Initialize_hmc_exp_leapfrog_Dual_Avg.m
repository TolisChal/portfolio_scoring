function [eps_step, X0, prob_step_in] = Initialize_hmc_exp_leapfrog_Dual_Avg(A, b, X0, c, T, N, delta)
% Reflective Hamiltonian Monte Carlo to sample from f truncated to the set
% of the correlation matrices
    
    d = length(c);
    
    cc = c / T;

    Delta_max = 1000;
    
    eps_step = 0.01;
    mu = log(10*eps_step);
    log_tilde_eps = 0;
    H_tilde = 0;
    gamma = 0.05;
    t0 = 10;
    k = 0.75;
    
    total_num_steps = 0;
    total_num_steps_in = 0;
    
    %h = waitbar(0,'Estimating a good step length for ReHMC...');
    for i = 1:N
        %i
        %for jj = 1:W
            x_counting_total = 0;
            v = randn(d, 1);
            v_pl = v;
            v_min = -v;
            X_pl = X0;
            X_min = X0;
            %Ti_pl = T0i;
            %Ti_min = T0i;
            X_nx = X0;
            %Ti_nx = T0i;
            %v_nx = v;
            grad_x_pl = cc; 
            grad_x_min = grad_x_pl;
            %h1 = a * (X0'*sigma*X0 - q*(muu'*X0)) + 0.5 * (v' * v);
            h1 = cc' * X0 + 0.5 * (v' * v);
            %exp(-h1)
            uu = log(rand) - h1;
            %x_pl = [X0(lower); T0i];
            %x_min = [X0(lower); T0i];
            %u = rand * exp(-h1);
            j = -1;
            s = 1;
            alpha = 0;
            while(s == 1)
                j = j + 1;
                na = 2^(j+1)-1;
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
                    v = v - (eps_step/2) * grad_x;
                    T = eps_step;
                    
                    lambdas = (b - A*X) ./ (A*v);
                    [~, pos_max] = min(1./lambdas);
                    total_num_steps = total_num_steps + 1;
                
                    while(true)
                    
                        %grad_x
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
                        X = X + lambda* v;
                        %Ti = Ti + lambda * v((n+1):end);
                    
                        %reflevt the ray
                        v = v - (2*A(pos_max, :)*v)*A(pos_max, :)';
                        %v = v*norm_v;
                        %p = pos(pos_max);
                        %    vv(p) = -vv(p);
                        %    v(1:n) = vv;
                    
                        T = T - lambda;
                    end
                    grad_x = cc; 
                    v = v - (eps_step/2) * grad_x;
  
                    %hj = a * (X'*sigma*X - q*(muu'*X)) + 0.5 * (v' * v);
                    hj = cc' * X + 0.5 * (v' * v);
                    alpha = alpha + min(1, exp(-hj + h1));
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
                    X_nx = X_rnd_j;
                    %Ti_nx = Ti_rnd_j;
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
                
            end
            
            H_tilde = (1 - 1/(i+t0))*H_tilde + (1/(i+t0)) * (delta - alpha/na);
            log_eps = mu - (sqrt(i) / gamma) * H_tilde;
            log_tilde_eps = (mu^(-k))*log_eps + (1 - mu^(-k))*log_tilde_eps;
            eps_step = exp(log_eps);
            %h2 = f(X_nx, Ti_nx, f_utils) + 0.5 * (v_nx' * v_nx);
            %u = log(rand);
            %filter_num = filter_num + 1;
            %h1
            %h2
            %if (u < h1-h2)
                %filter_acc = filter_acc + 1;
            X0 = X_nx;
                %T0i = Ti_nx;
            %end
        %end
        
        %correlation_matrices{i} = X0;
        %partial_variances{i} = diag(T0i);
        %i
        %filter_acc
        %filter_num
        
        %waitbar(i/N);
    end
    %acc_prob = filter_acc / filter_num;
    eps_step = exp(log_eps);
    prob_step_in = total_num_steps_in / total_num_steps;
    %exp(log_tilde_eps)
    %close(h);
end

