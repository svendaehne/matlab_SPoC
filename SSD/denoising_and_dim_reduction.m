%% INTRO
% Example script that shows how to use SSD for denoising and 
% dimensionality reduction. 
% Please also read (and cite) the following paper for more information on
% dim-reduction and denoising with SSD:  
%
% Haufe, S., Dahne, S., & Nikulin, V. V. Dimensionality reduction for the 
% analysis of brain oscillations. NeuroImage, 101:583-597, 2014 
% DOI: 10.1016/j.neuroimage.2014.06.073



%% THE DATA
% We assume X contains the continous EEG/MEG data with 
% columns of X corresponding to channels and rows to time points, 
% i.e. size(X) = [n_samples, n_channels]. 
% If you want to apply cSPoC, X is a cell array with X{n} being
% the n'th dataset. The preprocessing has to be applied to each 
% dataset individually

X = randn(10000, 50); % no real data in this script, only random numbers

% parameters for SSD, see the doc-string of ssd.m for more information
bands = [10,12; 8,14;  9,13];
sampling_freq = 200;

% some random target variable for SPoC
z = randn(1, size(X,1)/sampling_freq);

%% SSD FOR DIMENSIONALITY REDUCTION
% Here we are performing explicit dimensionality reduction via SSD. That
% means we apply SSD first and then carry out additional optimization (e.g.
% SPoC) on the time courses of a subset of the SSD components. 
% In order to apply the SPoC filters on the original input data, we will
% map the SPoC spatial filters from the SSD space to the original input space. 

% First we call the SSD function which returns the 
% original data bandpassed and projected onto the SSD components.
% Also, in order to compute the spatial patterns of extracted 
% (m/c)SPoC sources,  we need the covariance matrix of the 
% bandpassed sensorspace data.
% See the SSD code for more information on the paramters

[W_ssd, A_ssd, lambda_ssd, C, X_ssd]  = ssd(X, bands, sampling_freq, [],[]); 


% Dimensionality reduction by choosing only the first few 
% SSD components

n_ssd_components = 15; 
X_ssd = X_ssd(:,1:n_ssd_components);
W_ssd = W_ssd(:,1:n_ssd_components);


% X_ssd now contains the time-courses of the first n ssd components.
% Next we apply any of the (m/c)SPoC methods (depends on what you want to do)
% in the dimensionality-reduced space. 
X_ssd_epoched = permute(reshape(X_ssd, [sampling_freq, length(z), n_ssd_components]), [1,3,2]);
W_spoc = spoc(X_ssd_epoched, z);


% Note that the spatial filters are now defined in SSD space and have 
% to be mapped to sensor space to be applicable to the original data.

W = W_ssd * W_spoc;


% Now the source time courses can be extracted from the data using
% the combined SSD and SPoC filters

S = X * W;


% The spatial patterns are computed using the covariance matrix of the
% bandpass filtered sensor space data

A =  C * W;


%% SSD FOR DENOISING 
% In this alternavtive way of using SSD, the dimensionality reduction is
% carried out implicitly. Like in the example above, we begin by computing
% the SSD filters and patterns. Then the data is projected into SSD space.
% After that, the data is projected back into the original input space
% (let's call it "channel space") by using only a subset of the SSD
% components in the backprojection. The backprojected data has the same
% shape (i.e. dimensionality) as the original data but, by leaving out some
% SSD components during the backprojection, it has lower rank and
% (hopefully) a higher SNR in the frequency band of interest. 


% First we call the SSD function. This time we are only interested in the
% patterns and filters, i.e. A and W (and maybe also the lambda values as a
% basis for deciding how many components to keep)

[W_ssd, A_ssd, lambda_ssd] = ssd(X, bands, sampling_freq, [], []);

% Using the denoising function, we first project the data onto all SSD
% components and then backproject only the first n SSD components.

n_ssd_components = 15;
X_bp = denoising_by_low_rank_factorization(W_ssd, A_ssd, X, n_ssd_components);

% Now we can continue with the backprojected data (X_bp) and, for example,
% apply SPoC on it. 
% But don't forget that (m/c)SPoC require band-pass filtered data. So X_bp
% must be band-pass filtered first!

[b,a] = butter(5, [8,12]/(sampling_freq/2)); 
X_bp_flt = filter(b,a, X_bp); % This is just an example. You can use your
                            % favorite FIR or IIR filter settings /
                            % function.
X_bp_flt_epoched = permute(reshape(X_bp_flt, [sampling_freq, length(z), size(X,2)]), [1,3,2]);

[W, A, lambda] = spoc(X_bp_flt_epoched, z);

% Now we have spatial filters W and spatial patterns A that are defined in
% channel space. But note that the filters are optimized to work on
% denoised data! This means that if you want to apply them to test data,
% you have to do the denoising with the previously computed SSD filters and
% patterns!


                            