function [Copulas_feats, Copulas_sum_feats]=comp_indicators_feats(Copulas)

[redu, redl, blueu, bluel] = indicator_mask_corners(Copulas{1});

Copulas_sum_feats=zeros(size(Copulas,2),4);
Copulas_feats=zeros(size(Copulas,2),nchoosek(4,2));

for i = 1:size(Copulas,2)
    Copulas_sum_feats(i,:)=[sum(sum(redl.*Copulas{i})) sum(sum(redu.*Copulas{i})) sum(sum(bluel.*Copulas{i})) sum(sum(blueu.*Copulas{i}))];
    c=combnk(Copulas_sum_feats(i,:),2);
    Copulas_feats(i,:)=(c(:,1)./c(:,2))';
end