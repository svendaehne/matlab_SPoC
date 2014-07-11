function [W, A, lambda, C_s, X_ssd] = ssd(X, freq, sampling_freq, filter_order, epoch_indices)
% SSD - Spatio-Spectral Decomposition
%
% [W, A, lambda, C_s, X_ssd] = ssd(X, freq, sampling_freq, filter_order, epoch_indices)
%
% This is a function for the extraction of neuronal oscillations 
% with optimized signal-to-noise ratio. The algorithm maximizes 
% the power at the center frequency (signal of interest) while simultaneously suppressing it
% at the flanking frequency bins (considered noise). 
% 
% INPUT: 
%     X -     a matrix of size TxM with T samples and M channels of 
%             raw EEG/MEG/LFP data
%     freq -  3 x 2 matrix with the cut-off frequencies. 
%             First row: cut-off frequencies for band-pass of the to be extracted 
%             oscillations.
%             Second row: cut-off frequencies for the slowest and highest 
%             frequencies defining flanking intervals.
%             Third row: cut-off frequencies for the band-stop filtering of 
%             the central frequency process.
%     sampling_freq -     the sampling frequency (in Hz) of X (the data)
%     filter_order  -     filter order used for butterworth bandpass and 
%                         bandstop filtering. If unsure about it, use [] then
%                         the default value of order = 2 will be used
%     epoch_indices -     a matrix of size N x 2, where N is the number of 
%                         good (i.e. artifact free) data segments. Each row of 
%                         this matrix contains the first and the last sample index of
%                         a data segment, i.e. epoch_indices(n,:) = [1000, 5000]
%                         means that the n'th segment starts at sample 1000 and
%                         ends at sample 5000. If all data is useable or if unsure,
%                         use [], then the default of [1, size(X,1)] will be used. 
%                         
% OUTPUT:
%     W -     the de-mixing matrix. Each column is a spatial filter and the
%             timecourse of the SSD components is extracted with X * W
%     A -     the spatial patterns (also called mixing matrix) with the i'th column
%             corrresponding to the i'th SSD component
%     lambda - the eigenvalues corresponding to each component. The stronger
%              the eigenvalue the better is the ratio between the signal and noise. 
%              The components are sorted in the descending order (first components 
%              have the largest SNR)
%      C_s -  the covariance matrix of X after bandpass filtering with the band
%             defined in freq(1,:)
%      X_ssd - the bandpass filtered data projected onto the SSD components, 
%              i.e. X_ssd = X_s * W, where X_s is the bandpass filtered version of X
%              
% 
% EXAMPLE:
%   Let us consider that we want to extract oscillations in the 10-12 Hz
% frequency range with sampling frequency 200 Hz, then the FREQ can be defined as
% FREQ=[10 12; 8 14; 9 13]. Here we want to extract oscillations in 10-12
% Hz range, and flanking noise is defined as band-pass filtered data in
% 8-14 Hz with the following band-stop filtering in 9-13 Hz in order 
% to prevent spectral leakage to flanking noise from the signal of interest
% (10-12 Hz in this case). ORDER=2 ( 4 with filtfilt) 
% We want only data from 2 seconds to 100 seconds and then from 110 seconds to 150 seconds
% then EPOCH=[2*SF 100*SF; 110*SF 150*SF]
% where SF=200. 
%
% The whole command is then written as:
% [W, A, lambda, C_s, X_ssd] = ssd(X, FREQ, SF, ORDER, EPOCH); 
%
%
% References:
%
% Nikulin VV, Nolte G, Curio G. A novel method for reliable and fast extraction
% of neuronal EEG/MEG oscillations on the basis of spatio-spectral decomposition.
% NeuroImage, 2011, 55: 1528-1535.
%
% Haufe, S., Dahne, S., & Nikulin, V. V. Dimensionality reduction for the 
% analysis of brain oscillations. NeuroImage, 2014 (accepted)
%


%% check input arguments

% make sure FREQS has the correct dimensions
if not( size(freq,1)==3 && size(freq,2)==2 )
  error('FREQS must be a 3 by 2 matrix, i.e. three bands must be specified!');
end

% check that the given frequency bands are nested
signal_band = freq(1,:); % signal bandpass band
noise_bp_band = freq(2,:); % noise bandpass band
noise_bs_band = freq(3,:); % noise bandstop band
if not( noise_bs_band(1) < signal_band(1) && ...
        noise_bp_band(1) < noise_bs_band(1) && ...
        signal_band(2) < noise_bs_band(2) && ...
        noise_bs_band(2) < noise_bp_band(2) )
  error('The bands must be nested, i.e. the first band within the second and the second within the third!');
end

% default values for optional arguments
if isempty(filter_order)
    filter_order = 2;
end
if isempty(epoch_indices)
    epoch_indices = [1, size(X,1)];
end

% indices of good segments
ind=[];
for n=1:size(epoch_indices,1)
    ind=[ind epoch_indices(n,1):epoch_indices(n,2)];
end

%% filtering of data

% Creating filters
[b,a]=butter(filter_order, signal_band/(sampling_freq/2));
[b_f,a_f]=butter(filter_order, noise_bp_band/(sampling_freq/2));
[b_s,a_s]=butter(filter_order, noise_bs_band/(sampling_freq/2),'stop');


% Covariance matrix for the center frequencies (signal)
X_s = filtfilt(b,a,X);
C_s = cov(X_s(ind,:),1);

% Covariance matrix for the flanking frequencies (noise)
X_tmp = filtfilt(b_f,a_f,X);
X_tmp = filtfilt(b_s,a_s,X_tmp);
C_n = cov(X_tmp(ind,:),1);
clear X_tmp

%% Generalized eigenvalue decomposition

[W,D]= eig(C_n,C_s+C_n);
D=1-diag(D);
A=inv(W)'; % A is a matrix with the patterns (in columns)

[lambda,sorted_idx]=sort(D,'descend'); % sorting eigenvalues
W=W(:,sorted_idx);
A=A(:,sorted_idx);

%% apply SSD filters to the data

if nargout > 4
    X_ssd = X_s * W;
end



