function  [blue_mask, red_mask]=indicator_mask(X, band_ratio)

if nargin==1
    band_ratio=0.2;
end

A=ones(size(X));

band=round(band_ratio*size(X,1));

M=((tril(A,-band)==0+triu(A,band)==0)==0);

blue_mask=M-(M+flip(M)==2);

red_mask=flip(M)-(M+flip(M)==2);
