function X_denoised = denoising_by_low_rank_factorization(W, A, X, n_components)
% 
% Denoising by low-rank factorization. Computes a low-rank approximation of
% the input data, given a set of spatial filters and corresponding spatial
% patterns. 
% 
% X_denoised = denoising_by_low_rank_factorization(W, A, X, n_components)
%
% The input data X are projected onto a subset of components using the
% filter matrix W and then this subset is backprojected using the patterns
% in A. By using only a subset of components, the data is essentially
% "denoised" without changing the dimensionality. The rank is affected
% however!
%
% Input:
% W - matrix of spatial filters (for example obtained from SSD or ICA, etc).
%         Each column is one filter.
% A - matrix of spatial patterns (for example obtained from SSD or ICA, etc).
%         Each column is one pattern.
% X - the input data, size(X) can either be [n_samples, n_channels] or 
%         [n_samples, n_channels, n_epochs]
% n_components - number of components to use from W and A
%
% Output:
% X_denoised - same format as X but containing only the first n components
%               of W and A
%
%
% See the section on low-rank factorization in Haufe et al. 2014 for more
% for more details on the approach.
%
% Ref:  
% Haufe, S., Dähne, S., & Nikulin, V. V. (2014). "Dimensionality 
% reduction for the analysis of brain oscillations". NeuroImage, 101, 583-597.
% 
% sven.daehne@tu-berlin.de, 2015



%% some checks on the input

if isempty(W) && isempty(A)
    error('Both W and A are empty! At least one of them has to given.')
end

if isscalar(n_components)
    idx = 1:n_components;
else
    idx = n_components;
end

% make sure X is in the correct format
x_was_epoched = false;
if ndims(X) == 3
    % we assume X is epoched, i.e. size(X) = [T_epo, N_channels, N_epo]
    % here we make it 'continous'
    [T_epo, Nx, Ne] = size(X);
    X = reshape( permute(X, [1,3,2]), [T_epo*Ne, Nx]);
    x_was_epoched = true;
end
Nx = size(X,2);

% if the patterns are not given, compute them using the filter/pattern
% transformation
if isempty(A)
    C = cov(X);
    A = C * W / (W' * C * W);
end

% if the filters are not given, we have to assume that X is already in
% component space, i.e. that the spatial filters have been applied to X 
% outside this function
if isempty(W)
    W = eye(Nx);
end
    
%% low-rank approximation 

W = W(:,idx);
A = A(:,idx);

X_denoised = X * W * A';

%% 
if x_was_epoched
    % bring X_denoised in the same epoch format that X was in
    X_denoised = permute(reshape(X_denoised, [T_epo, Ne, Nx]), [1,3,2]);
end