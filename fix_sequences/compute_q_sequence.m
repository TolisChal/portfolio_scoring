function [qs, ptfs, q_min_gl, q_max, vols] = compute_q_sequence(sigma, mu, min_vol, max_vol, N)

    n = length(mu);
    
    A = [ eye(n) ; -eye(n)];
    b = [ ones(n,1) ; zeros(n,1)];
    Aeq = [ ones(1,n) ];
    beq = [ 1 ];

    options = optimoptions('fmincon','MaxFunctionEvaluations', 20000, 'Display', 'off', 'ConstraintTolerance', 1e-10, 'MaxIterations',5000, 'OptimalityTolerance', 1e-10);
    
    q_min_gl = 0;
    q_min = 0;
    q_max = 5;
    q_max_prev = 0;
    Ptf_q = (1/n)*ones(n,1);
    
    vol_prev = 0;
    max_vol_gl = max_vol;
    
    while (true)
        
        Ptf = fmincon(@(x) target_q_volatility_fun(x, sigma, mu, q_max), Ptf_q, A, b, Aeq, beq, [], [], [], options);
        vol = Ptf'*sigma*Ptf;
        
        if (vol >= max_vol_gl)
            max_vol = vol;
            break;
        elseif (abs(vol - vol_prev)/vol < 0.001)
            q_max = q_max/1.5;
            q_max_prev = q_max_prev / 1.5;
            max_vol = vol;
            break;
        else
            q_max_prev = q_max;
            q_max = 1.5*q_max;
            vol_prev = vol;
        end
        
    end
    q_prev = 0;
    q_max_gl = q_max;
    q_max_prev_gl = q_max_prev;
    while (true)
        
        q = (q_max+q_max_prev)/2
        
        Ptf = fmincon(@(x) target_q_volatility_fun(x, sigma, mu, q), Ptf_q, A, b, Aeq, beq, [], [], [], options);
        vol = Ptf'*sigma*Ptf;
        %abs(vol - max_vol)/max_vol 
        if (abs(vol - max_vol)/max_vol <= 0.001 & abs(q-q_prev)/q < 0.001)
            q_max = q;
            max_vol = vol;
            break;
        elseif (vol < max_vol && abs(vol - max_vol)/max_vol > 0.001)
            q_max_prev = q;
        else
            q_max = q;
        end
        q_prev = q;
    end
    
    step = (max_vol - min_vol) / (N+1);
    vols = min_vol:step:max_vol
    %vols = vols(2:end-1);

    counter = 1;
    qs = zeros(1,N);
    ptfs = zeros(n, N);
    for i = vols
        
        q_min_iter = q_min;
        q_max_iter = q_max;
        while(true)
            
            q = (q_min_iter + q_max_iter) / 2;
        
            Ptf = fmincon(@(x) target_q_volatility_fun(x, sigma, mu, q), Ptf_q, A, b, Aeq, beq, [], [], [], options);
            vol = Ptf'*sigma*Ptf;
            
            if (abs(vol-i)/i < 0.001)
                qs(counter) = q;
                ptfs(:, counter) = Ptf;
                counter = counter + 1;
                vol
                break;
            elseif (vol < i)
                q_min_iter = q;
            else
                q_max_iter = q;
            end
            
        end
        q_min = q;
    end
    qs
    q_min_gl = qs(1);
    q_max = qs(end);
    
    qs = qs(2:(end-1));
    ptfs = ptfs(:,2:(end-1));
    
    ptfs
end