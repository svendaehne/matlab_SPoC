function spoc_example()
% A small example that illustrates how SPoC works and how it can be
% applied. In this example we assume that we have a certain number of
% channels (Nx)that measure a mixture of source activity. The epoch-wise
% variance of one of the sources (target source) is the target variable z.
%
% Note that there is training and test data! 
%
% Run this script a few times to get a feeling for the influence of the
% paramters (SNR, Ns, Ne_tr, etc).


%% parameters

SNR = 0.2; % signal-to-noise ratio (between 0 and 1) in terms of variance explained by the target source
Ns = 5; % number of sources, the first one is the 'target source'
Nx = 15; % number of simulated EEG sensors, 
Ne_tr = 50; % number of training epochs/trials
Ne_te = 50; % number of test epochs/trials
Ne = Ne_tr + Ne_te;
Te = 200; % number of samples per epoch
samples_per_second = 100;

tr_idx = 1:Ne_tr;
te_idx = (1:Ne_te) + Ne_tr;

% make sure the SNR is between 0 and 1
SNR = max(0,SNR);
SNR = min(1,SNR);


%% data in source space

S = randn(Te*Ne, Ns);
[b,a] = butter(5, [3,4]/samples_per_second*2);
S = filtfilt(b, a, S);
S_env = abs(hilbert(S));

z = S_env(:,1).^2;  % Target function for SPoC analysis: the actual 
                    % squared envelope (i.e. power time course) of one of
                    % the sources. The task is to recover exactly that
                    % source.

% plot source signals
figure
plot_time_courses(my_zscore(S), 1,'s', 1);
title({'data in source space', 's_1 is the target source'})
xlabel('time')
set(gca, 'xtickLabel',[])
ylabel('sources')

%% mix the sources

% create some random smooth source patterns
A = randn(Nx,Ns);
A = filtfilt([1,1], 1, A);

% project sources to sensors (i.e. mix the sources)
X_s = S(:,1) * A(:,1)';
X_bg = S(:,2:end) * A(:,2:end)';

X_s = X_s ./ norm(X_s(:),'fro');
X_bg = X_bg ./ norm(X_bg(:),'fro');
X = SNR*X_s + (1-SNR)*X_bg;

% add some small sensor noise
X_noise = randn(size(X));
X_noise = X_noise/norm(X_noise(:),'fro');
X = X + 0.05*X_noise;

% plot sensor signals
figure
plot_time_courses(my_zscore(X), 1,'x');
title('data in sensor space')
xlabel('time')
set(gca, 'xtickLabel',[])
ylabel('sensors')


%% reshape sensor signals into epochs

X_epo = permute(reshape(X, [Te, Ne, Nx]), [1,3,2]);
z_epo = mean(reshape(z, [Te, Ne]));

% split training and test data
X_tr = X_epo(:,:,tr_idx);
z_tr = z_epo(tr_idx);

X_te = X_epo(:,:,te_idx);
z_te = z_epo(te_idx);

%% unmix with SPoC, only using training data

[W,A_est,lambda] = spoc(X_tr,z_tr);

%% project data to first SPoC component

s_est = zeros(Te, Ne);
for k=1:Ne
    s_est(:,k) = squeeze(X_epo(:,:,k)) * W(:,1);
end
p_est = var(s_est);


%% compare correlations obtained using channelwise power versus using SPoC component power

P = squeeze(var(X_te)); % channel-wise power (variance)
r = corrcoef([z_te',P']);
r = r(1,2:end);

r_spoc_tr = corrcoef(z_tr',p_est(tr_idx)');
r_spoc_tr = r_spoc_tr(1,2);
r_spoc_te = corrcoef(z_te',p_est(te_idx)');
r_spoc_te = r_spoc_te(1,2);

%% plot lambdas

figure
bar(lambda)
title('lambda values, i.e. SPoC optimization criterion')
ylabel('lambda')
xlabel('SPoC component')
xlim([0 length(lambda)+1])


%% plot time course and envelope of true and estimated source
figure
rows = 1;
cols = 5;
subplot(rows,cols,1:(cols-1))
s_true = my_zscore(S(:,1));
s_est_cnt = my_zscore(reshape(s_est, [Te*Ne, 1]));
s_est_cnt = s_est_cnt * sign(s_true'*s_est_cnt);

hold on
h = plot([s_true, s_est_cnt]);
% plot(abs(hilbert(([s_true, s_est_cnt]))))
legend(h, {'true source signal', '1st SPoC component'})
title('time course and envelope of true source and estimated source (1st SPoC component)')
ylabel('normalized (z-scored) amplitude [a.u.]')
xlabel('time')
set(gca, 'xtick',[]);
box on

% scatter true and estimated power
subplot(rows,cols,cols)
scatter(my_zscore(z_te), my_zscore(p_est(te_idx)))
xlim([-3,3])
ylim([-3,3])
box on
title({'true source power vs estimated source power', 'on test data'})
axis square
xlabel('power of 1st SPoC component')
ylabel('power of true source')

%% plot channelwise power correlations vs SPoC correlations

% correlations
figure
hold on
box on
bar(r, 'k')
plot(ones(size(r))*r_spoc_tr,'--r')
plot(ones(size(r))*r_spoc_te,'r')
title({'correlation of power with target function'})
xlim([0,Nx+1])
ylim([-0.2,1.3])
set(gca, 'ytick', -1:0.2:1)
xlabel('channel index')
ylabel('r')
legend({'corr z with channel-wise power', ...
    'corr z with 1st SPoC component power (training data)', 'corr z with 1st SPoC component power (test data)'}, ...
    'location','best')

