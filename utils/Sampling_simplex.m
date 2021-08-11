function [ E ] = Sampling_simplex( N, m, method )
%SAMPLING is sampling uniformly from a unit simplex
%
% N: number of dimension of the simplex (-1 ?)
% m: number of random vectors
% method: method used 'RM' = Rubinstein and Melamed, 'RK' = Rubinstein and Kroese
% Ref.: - Rubinstein and Melamed (1998),  Algorithm 2.7.1
%       - Rubinstein and Kroese (2007), Algorithm 2.5.3

    switch method
        case 'RM' % Rubinstein and Melamed
            pd = makedist('Exponential');
            Y = random(pd, N, m);
            T = sum(Y);
            E = Y./(ones(N,1)*T);
        case 'RK' % Rubinstein and Kroese
            Y = zeros(N+1,m);
            Y(N+1,:) = 1;
            Y(2:N,:) = sort(rand(N-1,m));
            E = Y(2:end, :)-Y(1:end-1, :);
    end

end

