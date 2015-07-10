function [X, Y, Sx, Sx_env, Sy, Ax, Ay] = create_mspoc_example_data(K, n_sources, ...
            n_channels, HRF, SNR, samples_per_second, Ne, Te, band)

if length(n_sources) == 1
    n_sources = [1,1]*n_sources;
end

N = 10;
Ne = Ne + N;

%% create the source signals

% oscillatory source signals and envelopes -> X-modality
Sx = randn(Te*Ne, n_sources(1));
[b,a] = butter(5, band/samples_per_second*2);
Sx = filtfilt(b, a, Sx);
Sx_env = abs(hilbert(Sx));
Sx = Sx ./ Sx_env;
[b,a] = butter(5, 0.25/samples_per_second*2);
Sx_env = filtfilt(b, a, Sx_env);
Sx = Sx .* Sx_env;
 
% slowly varying sources -> Y-modality
Sy = randn(Te*Ne, n_sources(2));
[b,a] = butter(5, band/samples_per_second*2);
Sy = filtfilt(b, a, Sy);
Sy = abs(hilbert(Sy));
[b,a] = butter(5, 0.25/samples_per_second*2);
Sy = filtfilt(b, a, Sy);

% couple the first K sources in the two modalities by setting the envelopes
% of the first k X-sources to be the time-courses of the first k Y-sources
Sx(:,1:K) = Sx(:,1:K) ./ Sx_env(:,1:K);
Sx(:,1:K) = Sx(:,1:K) .* Sy(:,1:K);
Sx_env = abs(hilbert(Sx));


% implements sub-sampling by averaging
Sy = reshape(Sy, [Te, Ne, n_sources(2)]);   
Sy = squeeze(mean(Sy));

% filter the Y-source time-courses with the HRF
Sy = filter(HRF, 1, Sy);

Sx = Sx((N*Te+1):end,:);
Sx_env = Sx_env((N*Te+1):end,:);
Sy = Sy((N+1):end,:);


%% mix the sources

S = {Sx, Sy};
X = cell(1,2);
A = cell(1,2);
for n=1:2
    % create some random smooth source patterns
    A{n} = randn(n_channels(n),n_sources(n));
%     d = ceil(n_channels(n)*0.01);
%     b = ones(1,max(d,2));
    b = ones(1,3);
    A{n} = filtfilt(b/sum(b), 1, A{n});
    
    % project sources to sensors (i.e. mix the sources)
    X_s = S{n}(:,1:K) * A{n}(:,1:K)';
    X_bg = S{n}(:,(K+1):end) * A{n}(:,(K+1):end)';
    
    X_s = X_s ./ norm(X_s(:),'fro');
    X_bg = X_bg ./ norm(X_bg(:),'fro');
    X{n} = SNR*X_s + (1-SNR)*X_bg;
    
    % add some small sensor noise
    X_noise = randn(size(X{n}));
    X_noise = X_noise/norm(X_noise(:),'fro');
    X{n} = X{n} + 0.005*X_noise;
end

Y = X{2};
X = X{1};
Ax = A{1};
Ay = A{2};
