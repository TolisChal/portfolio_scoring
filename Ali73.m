function [Score] = Ali73(AssetRet,Ret)
%ALI73 Summary of this function goes here
%   Detailed explanation goes here
    
    U = AssetRet - Ret;
    
    Y = U(U>=0);
    K = length(Y);
    
    X = U(U<0);
    J = length(X);
    
    A = zeros(length(Y)+1,1);
    A(1) = 1;
    
    for h=1:J
        for k = 1:K
            A(k+1) = ( Y(k)*A(k+1)-X(h)*A(k) )/( Y(k)-X(h));
        end 
    end
    
    Score = A(end);
    
end

