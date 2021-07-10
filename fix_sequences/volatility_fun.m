function res = volatility_fun(Ptf,sigma)
    res = Ptf'*sigma*Ptf;
end

