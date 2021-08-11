function y = low_risk(x)
    
    x = 0.5 +(0.5 - x);
    y = 1.031*exp(x)-0.8;

end