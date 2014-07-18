function cspoc_example()
% A small example that illustrates how cSPoC works and how it can be
% applied.
% In this example we assume that we have two data sets X and Y, each with
% a certain number of channels (Nx, Ny) that measure a mixture of source 
% activity. The envelope of a source in X is coupled to a source in Y,
% these are the "target sources" that we seek to extract with cSPoC.
%
% Note that there is training and test data! 
%
% Run this script a few times to get a feeling for the influence of the
% paramters (SNR, Ns, Ne_tr, etc).


%% parameters

N = 3; % number of datasets
K = 2; % number of "target source pairs", i.e. sources with coupled envelopes
n_sources = K + 1 + (1:N); % number of total simulated sources per dataset
n_channels = max(5,K) + (1:N); % number of channels per dataset

SNR = 0.2; % signal-to-noise ratio (between 0 and 1) in terms of variance explained by the target source

Ne_tr = 50; % number of training epochs/trials
Ne_te = 50; % number of test epochs/trials
Ne = Ne_tr + Ne_te;
Te = 100; % number of samples per epoch
samples_per_second = 100;

tr_idx = 1:Ne_tr;
te_idx = (1:Ne_te) + Ne_tr;

% make sure the SNR is between 0 and 1
SNR = max(0,SNR);
SNR = min(1,SNR);

S = cell(1,N); % source signals (time courses)
S_env = cell(1,N); % envelopes of source signals
X = cell(1,N); % sensor signals


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


%% data in source space

fprintf('creating data \n')

% generate some random oscillations
for n=1:N
    S{n} = randn(Te*Ne, n_sources(n));
    [b,a] = butter(5, [3,4]/samples_per_second*2);
    S{n} = filtfilt(b, a, S{n});
    S_env{n} = abs(hilbert(S{n}));
end

% make K source envelopes correlate
for n=2:N
    S{n}(:,1:K) = S{n}(:,1:K) ./ S_env{n}(:,1:K);
    S{n}(:,1:K) = S{n}(:,1:K) .* S_env{1}(:,1:K);
end
% plot source signals
figure
for n=1:N
    subplot(N,1,n)
    plot_time_courses(zscore(S{n}), 1,'s',1:K);
    title(sprintf('dataset #%d in source space',n))
    xlabel('time')
    set(gca, 'xtickLabel',[])
    ylabel('sources')
end
pause(0.1)

%% mix the sources

A = cell(1,N);
for n=1:N
    % create some random smooth source patterns
    A{n} = randn(n_channels(n),n_sources(n));
    A{n} = filtfilt([1,1], 1, A{n});
    
    % project sources to sensors (i.e. mix the sources)
    X_s = S{n}(:,1:K) * A{n}(:,1:K)';
    X_bg = S{n}(:,(K+1):end) * A{n}(:,(K+1):end)';
    
    X_s = X_s ./ norm(X_s(:),'fro');
    X_bg = X_bg ./ norm(X_bg(:),'fro');
    X{n} = SNR*X_s + (1-SNR)*X_bg;
    
    % add some small sensor noise
    X_noise = randn(size(X{n}));
    X_noise = X_noise/norm(X_noise(:),'fro');
    X{n} = X{n} + 0.01*X_noise;
end


% plot sensor signals
figure
for n=1:N
    subplot(N,1,n)
    plot_time_courses(zscore(X{n}), 1,'x');
    title(sprintf('dataset #%d in sensor space', n))
    xlabel('time')
    set(gca, 'xtickLabel',[])
    ylabel('sensors')
end
pause(0.5)

%% reshape sensor signals into epochs

X_epo = cell(1,N);
X_epo_tr = cell(1,N);
X_epo_te = cell(1,N);

for n=1:N
    % create epoched data structure
    X_epo{n} = permute(reshape(X{n}, [Te, Ne, n_channels(n)]), [1,3,2]);
    
    % split training and test data
    X_epo_tr{n} = X_epo{n}(:,:,tr_idx);
    X_epo_te{n} = X_epo{n}(:,:,te_idx);
end

%% unmix with cSPoC, only using training data
fprintf('\n\nbegin extracting envelope-correlated sources with cSPoC ... \n')

[W, A] = cspoc(X_epo_tr,1, ...
    'n_component_sets', K+1, ... % number of component pairs to be extracted
    'average_over_epochs', 0, ...
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
    ylim([-4,5])
    box on
    legend(lh,legend_string)
    title(sprintf('time-course of extracted source pair #%d', k))
    
    
    plot(Ne_tr*Te*[1,1], ylim, 'color',0.6*[1,1,1])
    text(Ne_tr*Te*1.01, 4, '--> test data')
    h=text(Ne_tr*Te*0.99, 4, 'training data <--');
    set(h, 'HorizontalAlignment', 'right');
    
end    
