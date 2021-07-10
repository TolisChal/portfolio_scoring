function res = neg_volatility_fun(Ptf,sigma)
    res = - Ptf'*sigma*Ptf;
end

