function [Score] = Ali73_vecRet(AssetRet,Ptf)
%ALI73 Summary of this function goes here
%   Detailed explanation goes here
    
    NbObs       = size(AssetRet, 1);
    NbAssets    = size(AssetRet ,2);
    Score       = nan(NbObs,1);
    Vec_Obs     = 1:NbObs;
    
    U           = (AssetRet - AssetRet*Ptf*ones(1, NbAssets))';

    vec_K       = sum(U>=0);
    vec_J       = sum(U<0);
    %uniq_KJ     = unique([vec_K' vec_J'], 'rows');
    %nb_uniq_KJ  = size(uniq_KJ, 1);
    
    
    %for i = 1:nb_uniq_KJ
        
        %K = uniq_KJ(i,1);
        %J = uniq_KJ(i,2);
    for i=0:NbAssets
        
        K = i;
        J = NbAssets-i;
        BinPresence = (vec_K==K & vec_J==J);
        NbSelObs = sum(BinPresence);
        if(NbSelObs>0)
        
            SelObs      = Vec_Obs(vec_K==K & vec_J==J);
            NbSelObs    = length(SelObs);

            U_Sel       = U(:,SelObs);

            Y_tmp       = U_Sel( U(:,SelObs)>=0);    
            Y           = reshape(Y_tmp, [K NbSelObs])';

            X_tmp       = U_Sel(U(:,SelObs)<0);    
            X           = reshape(X_tmp, [J NbSelObs])';

            A       = zeros(NbSelObs, K+1);
            A(:,1)  = 1;

            for h=1:J
                for k = 1:K
                    A(:,k+1) = ( Y(:,k).*A(:,k+1)-X(:,h).*A(:,k) )./( Y(:,k)-X(:,h));
                end 
            end     

            Score(SelObs) = A(:,end);
            clear K J SelObs NbSelObs U_Sel Y_tmp Y X_tmp X A;
        end
    end
end

