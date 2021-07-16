function y = medium_risk(x)

    y = 8*(exp(-(x - 0.5).^2)-exp(-0.25)) + 0.2;

end