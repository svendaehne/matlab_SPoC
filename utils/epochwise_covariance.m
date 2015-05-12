function Cxxe = epochwise_covariance(X_epo)
%
% Cxxe = epochwise_covariance(X_epo)
%
% Computes epochwise covariance matrices.
%
% IN:
%   X_epo - epoched data with size(X_epo) = [T, n_channels, n_epochs]
%
% OUT:
%   Cxxe - timeseris of covariance matrices with 
%          size(Cxxe) = [n_channels, n_channels, n_epochs]

% 2014, sven.daehne@tu-berlin.de

n_channels = size(X_epo,2);
n_epos = size(X_epo,3);

Cxxe = zeros(n_channels, n_channels, n_epos);
for n=1:n_epos
    Cxxe(:,:,n) = cov(X_epo(:,:,n));
end