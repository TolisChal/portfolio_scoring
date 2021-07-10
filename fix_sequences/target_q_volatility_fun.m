function res = target_q_volatility_fun(Ptf,sigma,mu,q)
    res = Ptf'*sigma*Ptf - q*mu'*Ptf;
end

