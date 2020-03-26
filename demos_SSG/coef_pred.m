function T = coef_pred(train_lab)
% 2013-3-19
nClass = length(unique(train_lab));
train_size = length(train_lab);
T = zeros(nClass, train_size);
cls = unique(train_lab);
for k = 1 : nClass,
    c = cls(k);
    T(k, train_lab==c) = 1;
    %T(k, train_lab==c) = 1/sum(train_lab==c);
end
end