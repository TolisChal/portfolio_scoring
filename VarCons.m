function [c,ceq] = VarCons(Ptf,sigma, VarCst)
    ceq = Ptf'*sigma*Ptf - VarCst;
    c = [];
end

