function indicators = comp_indicators_mask(C)

[blue_mask, red_mask]=indicator_mask(C{1});

indicators=zeros(size(C,2),1);

for i = 1:size(indicators,1)
    indicators(i)=sum(sum(red_mask.*C{i}))/sum(sum(blue_mask.*C{i}));
end