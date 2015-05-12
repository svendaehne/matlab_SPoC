function rm_idx = detect_max_eigenvalue_epochs(X_epo, threshold, show_plot,perc_factor)
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

% detect outliers
if not(exist('threshold','var')) || isempty(threshold)
    p = percentiles(ev, [50,95]);
    if not(exist('perc_factor','var'))
        perc_factor = 2.5;
    end
    threshold = p(1) + perc_factor*(p(2)-p(1));
end

rm_idx = find(ev > threshold);



if not(exist('show_plot','var'))
    show_plot = true;
end

% plot the results
if show_plot
    figure
    
    hold on
    plot(ev)
    plot(1:Ne,threshold, '--r')
    if not(isempty(rm_idx))
        plot(rm_idx, threshold, '*k')
    end
    xlim([1,Ne])
    box on
    xlabel('epochs')
end




