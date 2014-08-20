function px = get_var_features(W, Cxxe)
% Computes the power timeseries (a.k.a. variance features) for each of the
% spatial filters contained in W, given the trialwise covariance matrix
% Cxxe. 
%
%

[Nx, ~, Ne] = size(Cxxe);
% vectorize the trialwise covariance time series
Cxxe_vec = reshape(Cxxe, [Nx*Nx, Ne]);
% vectorize the outer product of each column in W
W_vec = zeros(Nx*Nx, size(W,2));
for k=1:size(W,2)
    W_vec(:,k) = reshape(W(:,k)*W(:,k)', [Nx*Nx, 1]);
end
% power of projected x signal
px = W_vec' * Cxxe_vec;
