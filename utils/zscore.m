function z = zscore(x)
%ZSCORE Normalizes the columns of x to have zero mean and unit variance. 
z = (x-repmat(mean(x), [size(x,1),1])) ./ repmat(std(x), [size(x,1),1]);

