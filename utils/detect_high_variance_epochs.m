function rm_idx = detect_high_variance_epochs(X_epo, threshold, show_plot, perc_factor)
% 
% rm_idx = detect_high_variance_epochs(X_epo, threshold, show_plot)
%
% Displays epoch-wise variance and detects epochs with variance higher than
% a given or internally computed threshold. This is useful for the
% detection of artifact epochs and their removal. 
%
% In:
%   X_epo - data matrix of size(X_epo) = [T, n_channels, n_epos]
%           This input can also be two dimensional. In this case, it is
%           assumed that size(X_epo) = [n_channels, n_epos]
%
% Optional input:
%   threshold - threshold for mean channel-wise power, default: computed
%               based on percentiles of the channelwise averaged power
%   show_plot - show variance plot or not, default: 1 (yes)
%
% Out:
%   rm_idx - indices of epochs that exceed the threshold
%
%
% Examples:
%
% rm_idx = detect_high_variance_epochs(X_epo);
%
% rm_idx = detect_high_variance_epochs(X_epo, threshold, show_plot);
%
% rm_idx = detect_high_variance_epochs(squeeze(log(var(X_epo))), threshold, show_plot)



% 2015, sven.daehne@tu-berlin.de


if ndims(X_epo) == 3
    V = squeeze(var(X_epo));
else
    V = X_epo;
end

if not(exist('threshold','var')) || isempty(threshold)
    p = percentiles(mean(V), [50,95]);
    if not(exist('perc_factor','var'))
        perc_factor = 2.5;
    end
    threshold = p(1) + perc_factor*(p(2)-p(1));
end

if not(exist('show_plot','var')) || isempty(show_plot)
    show_plot = true;
end

rm_idx = find(mean(V,1) > threshold);


% plot the results
if show_plot
    rows = 4;
    cols = 6;
    figure
    
    tmp = repmat(1:(cols-1), [rows-1,1]) + repmat((0:(rows-2))'*cols, [1,cols-1]);
    subplot(rows,cols,tmp(:));
    imagesc(V)
    axis xy
    ylabel('channels')
    
    subplot(rows,cols, (1:(cols-1))+(rows-1)*cols )
    hold on
    plot(mean(V,1))
    plot(1:size(V,2),threshold, '--r')
    if not(isempty(rm_idx))
        plot(rm_idx, threshold, '*k')
    end
    xlim([1,size(V,2)])
    box on
    xlabel('epochs')
    
    subplot(rows,cols, cols:cols:((rows-1)*cols) )
    plot(mean(V,2),1:size(V,1))
    ylim([1,size(V,1)+0.0001])
    box on
    
end




