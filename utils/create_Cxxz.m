function Cxxz = create_Cxxz(Cxxe, z)
% Creates the z-weighted covariance matrix from a time series of covariance
% matrices given in Cxxe 

[Nx, ~, Ne] = size(Cxxe);
% vectorize the trialwise covariance time series
Cxxe_vec = reshape(Cxxe, [Nx*Nx, Ne]);
% compute the z-weighted average covariance matrix via matrix-vector multiplication
Cxxz = reshape(Cxxe_vec * z', [Nx, Nx]) / Ne;
 