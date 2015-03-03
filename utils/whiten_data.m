
function [X_white, M, M_inv] = whiten_data(X, min_var_explained, C)
% Returns the whitened data, the whiteneing matrix and its inverse
%
% X -   input data; is assumed to have dimensions: [samples x channels x epochs]
% MIN_VAR_EXPLAINED -   Between 0 and 1. If smaller than 1, dim reduction will be done.
% C -   covariance matrix of the data. Can be empty, then it will be computed
%       from X.

if not(isempty(X))
    X_orig = X;
    [T, n_channels, n_epos] = size(X);
    mu = squeeze(mean(mean(X,1),3));
end

if isempty(C)
    % reshape data to have dims [samples x channels] and remove the mean
    X = X(:,:,mask==1);
    X = reshape(permute(X, [1,3,2]), [T*sum(mask==1), n_channels]);
    X = X - repmat(mu, [T*sum(mask==1), 1]);
    C = cov(X);
end

% perform PCA and compute whitening/dim-reduction matrix
[V, D] = eig(C);
[ev_sorted, sort_idx] = sort(diag(D), 'descend');
V = V(:,sort_idx);
D = diag(ev_sorted);
% compute an estimate of the rank of the data
tol = ev_sorted(1) * 10^-6;
r = sum(ev_sorted > tol);
% compute cumulative variance explained
var_explained = cumsum(ev_sorted)/sum(ev_sorted);
n_components =  find(var_explained>=min_var_explained, 1);
n_components = min(n_components, r);
% construct the whitening matrix
M = diag(diag(D).^-0.5) * V'; % whitening filters are in the rows!!!
M_inv = V * diag(diag(D).^0.5); 
% dim-reduction
M = M(1:n_components, :);
M_inv = M_inv(:,1:n_components);

if not(isempty(X))
    % project the original data and do the dim-reduction
    X = X_orig;
    n_epos = size(X,3);
    X = reshape(permute(X, [1,3,2]), [T*n_epos, n_channels]);
    X = X - repmat(mu, [T*n_epos, 1]);
    X_white = (M * X')';
    X_white = permute(reshape(X_white, [T, n_epos, n_components]), [1,3,2]);
else
    X_white = [];
end