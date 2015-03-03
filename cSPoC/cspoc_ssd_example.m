%
% A small example that illustrates how cSPoC works and how it can be
% applied in conjunction with SSD (or any other linear method)
% preprocessing.
%
% In this example we assume that we have N datasets, each with
% a certain number of channels that measure a mixture of source 
% activity. The envelope of K sources in each dataset is coupled to a 
% corresponding sources in the other datasets,
% these are the "target sources" that we seek to extract with cSPoC.
%
% Note that there is training and test data! 
%
% Run this script a few times to get a feeling for the influence of the
% paramters (SNR, N, Ne_tr, etc).


%% parameters

N = 3; % number of datasets, e.g. subjects
K = 2; % number of "target source pairs", i.e. sources with coupled envelopes
n_sources = K + 10 + (1:N); % number of total simulated sources per dataset

n_channels = 100 + (1:N); % number of channels per dataset

n_ssd_components = 20; % number of SSD components used for dim-reduction

band = [8,12]; % frequency band of simulated data

SNR = 0.5; % signal-to-noise ratio (between 0 and 1) in terms of variance explained by the target source
Ne_tr = 250; % number of training epochs/trials
Ne_te = 250; % number of test epochs/trials
Ne = Ne_tr + Ne_te;
Te = 50; % number of samples per epoch
samples_per_second = 100;

tr_idx = 1:Ne_tr;
te_idx = (1:Ne_te) + Ne_tr;

% make sure the SNR is between 0 and 1
SNR = max(0,SNR);
SNR = min(1,SNR);

fprintf('\n')
fprintf('--- cSPoC example ---\n')
fprintf('\n')
fprintf('Number of simulated datasets = %d\n', N)
fprintf('Number of shared sources = %d\n', K)
fprintf('Number of total sources per dataset = ')
for n=1:(N-1); 
    fprintf(' %d, ', n_sources(n)); 
end
fprintf(' %d\n', n_sources(N));
fprintf('Number of channels per dataset = ')
for n=1:(N-1); 
    fprintf(' %d, ', n_channels(n)); 
end
fprintf(' %d\n', n_channels(N));

fprintf('\n')

%% create example data

fprintf('creating data \n')
[X, S, S_env, A_true] = create_cspoc_example_data(N, K, n_sources, n_channels,... 
                        SNR, samples_per_second, Ne, Te, band);

%% apply SSD for dimensionality reduction

fprintf('applying SSD for dimensionality reduction\n')
freq = [band; band+[-2,2]; band+[-1,1]];
for n=1:N

    % compute filters and patterns from SSD
    [W_ssd, A_ssd] = ssd(X{n}, freq, samples_per_second,[],[]);

    % denoise by projecting onto a subset of components
    X{n} = denoising_by_low_rank_factorization(...
            W_ssd, A_ssd, X{n}, n_ssd_components);
    
    % IMPORTANT!!!
    % If X is broadband, don't forget to band-pass filter it here! 
    %  X{n} = bandpass_filter(...);
    % In this example X is already band-limited because of the way it was
    % generated, so this step is not necessary.
    
end


%% reshape sensor signals into epochs

X_epo = cell(1,N);
X_epo_tr = cell(1,N);
X_epo_te = cell(1,N);

for n=1:N

    % create epoched data structure
    X_epo{n} = permute(reshape(X{n}, [Te, Ne, size(X{n},2)]), [1,3,2]);
    
    % split training and test data
    X_epo_tr{n} = X_epo{n}(:,:,tr_idx);
    X_epo_te{n} = X_epo{n}(:,:,te_idx);
end

%% unmix with cSPoC, only using training data

fprintf('\n\nbegin extracting envelope-correlated sources with cSPoC ... \n')

[W, A] = cspoc(X_epo_tr,1, ...
    'n_component_sets', K+1, ... % number of component pairs to be extracted
    'average_over_epochs', 1, ...
    'n_repeats', 5, ... % number of re-starts per component pair
    'verbose', 2 ...
    );

fprintf('\ndone!\n')


%% project data to cSPoC components and compute envelopes 

S_est = cell(1,N);
phi = cell(1,N);
for n=1:N
    S_est{n} = X{n} * W{n};
    phi{n} = abs(hilbert(S_est{n}));
end


%% plot extracted components

figure
rows = size(W{1},2);
cols = 1;

colors = colormap(lines(N));
lh = zeros(1,N);
legend_string = cell(1,N);
for k=1:rows
    subplot(rows,cols, k)
    hold on
    for n=1:N
        lh(n) = plot(S_est{n}(:,k), 'color', colors(n,:));
        plot(phi{n}(:,k), 'color', colors(n,:), 'linewidth',2)
        legend_string{n} = sprintf('source extracted from dataset #%d',n);
    end
    ylim([-5,6])
    box on
    legend(lh,legend_string)
    title(sprintf('time-course of extracted source pair #%d', k))
    
    
    plot(Ne_tr*Te*[1,1], ylim, 'color',0.6*[1,1,1])
    text(Ne_tr*Te*1.01, 4, '--> test data')
    h=text(Ne_tr*Te*0.99, 4, 'training data <--');
    set(h, 'HorizontalAlignment', 'right');
end    

%% plot spatial patterns
figure
rows = size(W{1},2);
cols = N;
for k=1:rows
    for n=1:N
        subplot(rows,cols, sub2ind([rows,cols],n,k))
        a_true = A_true{n}(:,k) / norm(A_true{n}(:,k));
        a = A{n}(:,k) / norm(A{n}(:,k));
        a = a * sign(corr(a, a_true));
        plot([a_true, a])
        axis tight
        xlabel('channel idx')
        ylabel('[a.u.]')
        title(sprintf('spatial pattern component %d, dataset %d', k, n))
        legend({'true pattern', 'estimated pattern'})
    end
end
