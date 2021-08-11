function res = FindInitPtf(Ptf,sigma, VarCst)
    res = abs(Ptf'*sigma*Ptf - VarCst);
end

