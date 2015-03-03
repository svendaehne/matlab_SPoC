function [X, S, S_env, A] = create_cspoc_example_data(N, K, n_sources, ...
            n_channels, SNR, samples_per_second, Ne, Te, band)


S = cell(1,N); % source signals (time courses)
S_env = cell(1,N); % envelopes of source signals
X = cell(1,N); % sensor signals


for n=1:N
    S{n} = randn(Te*Ne, n_sources(n));
    [b,a] = butter(5, band/samples_per_second*2);
    S{n} = filtfilt(b, a, S{n});
    S_env{n} = abs(hilbert(S{n}));
end

% make K source envelopes correlate
for n=2:N
    S{n}(:,1:K) = S{n}(:,1:K) ./ S_env{n}(:,1:K);
    S{n}(:,1:K) = S{n}(:,1:K) .* S_env{1}(:,1:K);
end

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
    X{n} = X{n} + 0.05*X_noise;
end
