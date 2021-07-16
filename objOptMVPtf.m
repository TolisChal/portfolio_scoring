function res = objOptMVPtf(Ptf, mu)

    meanPtf = mu'*Ptf;
    res = - meanPtf;
    %VolPtf  = sqrt(Ptf'*sigma*Ptf);
    %res = -meanPtf/VolPtf;
end

