function [X_epo, idx_flt, Cxxe] = median_filter_cov_timeseries(X_epo, N, show_plot)
% 
% rm_idx = detect_max_eigenvalue_epochs(X_epo, threshold, show_plot)
%
% Computes epochwise covariance matrices and their respective maximum
% eigenvalues. The function detects epochs with eigevalues higher than
% a given or internally computed threshold. This is useful for the
% detection of artifact epochs and their removal. 
%
% In:
%   X_epo - data matrix of size(X_epo) = [T, n_channels, n_epos]
%           If size(X_epo,1) == size(X_epo,2), it is assumed that X_epo is
%           already pre-computed time-series of covariance matrices 
%
% Optional input:
%   threshold - threshold for outlier detection, default: computed
%               based on percentiles of the eigenvalue distribution
%   show_plot - show variance plot or not, default: 1 (yes)
%
% Out:
%   rm_idx - indices of epochs that exceed the threshold
%
%
% Examples:
%
% rm_idx = detect_max_eigenvalue_epochs(X_epo);
%
% Cxxe = epochwise_covariance(X_epo);
% rm_idx = detect_max_eigenvalue_epochs(Cxxe);
%
% rm_idx = detect_max_eigenvalue_epochs(X_epo, threshold, show_plot);


% 2015, sven.daehne@tu-berlin.de


% epochwise covariance matrices
if size(X_epo,1) == size(X_epo,2)
    Cxxe = X_epo;
else
    Cxxe = epochwise_covariance(X_epo);
end

Ne = size(Cxxe,3);

% epochwise max eigenvalues 
ev = zeros(1,Ne);
for n=1:Ne
    ev(n) = max(eig(Cxxe(:,:,n)));
end

if not(exist('n','var')) || isempty(N)
    N = 3;
end

[ev_flt, idx_flt] = my_medfilt1(ev,N);

Cxxe = Cxxe(:,:,idx_flt);
X_epo = X_epo(:,:,idx_flt);

if not(exist('show_plot','var'))
    show_plot = true;
end

% plot the results
if show_plot
    figure
    
    hold on
    plot([ev', ev_flt])
    xlim([1,Ne])
    box on
    xlabel('epochs')
end




