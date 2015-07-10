

%% params


n_channels = [30,500];
K = 1; % number of "target source pairs", i.e. sources with coupled envelopes
n_sources = [15,15];
Nt = 0;

band = [8,12]; % frequency band of interest here
SNR = 0.7; % signal-to-noise ratio (between 0 and 1) in terms of variance explained by the target source

Ne_tr = 1000; % number of training epochs
Ne_te = 1000; % number of test epochs
Ne = Ne_tr + Ne_te;
Te = 100; % number of samples per epoch
samples_per_second = 200;

tr_idx = 1:Ne_tr;
te_idx = (1:Ne_te) + Ne_tr;

% make sure the SNR is between 0 and 1
SNR = max(0,SNR);
SNR = min(1,SNR);

if Nt == 0
    hrf = 1;
else
        hrf = exp(-(((0:Nt)-Nt/2).^2) / (2*(Nt/8)^2));
%     hrf_length_s = (Nt-1)*Te/samples_per_second;
%     hrf = spm_hrf(Te/samples_per_second, [0.15*hrf_length_s, 0.5*hrf_length_s, 1, 1, 6, 0, hrf_length_s])';
end

hrf = hrf/sum(hrf);

% figure
% plot(hrf)

%% create example toy data

[X, Y, Sx, Sx_env, Sy, Ax, Ay] = create_mspoc_example_data(K, n_sources, ...
            n_channels, hrf, SNR, samples_per_second, Ne, Te, band);

X = permute(reshape(X, [Te, Ne, n_channels(1)]), [1,3,2]);   
Y = Y';

Y_tr = Y(:,tr_idx);
X_tr = X(:,:,tr_idx);

%% optimize regularizers on training data

mspoc_params = struct('tau_vector',0:(length(hrf)-1)); 
[best_kappa_tau, best_kappa_y, xval_out] = ...
    optimize_mspoc_regularizers(X_tr, Y_tr, mspoc_params ...
        ,'kappa_tau_list', 0:0.2:1 ...
        ,'kappa_y_list', 0:0.2:1 ...
    );

%% mkSPoC on training data


% kappa_tau_list = 10.^(-4:0);
% kappa_y_list = 10.^(-4:0);

kappa_tau = 0; 10.^-2;
kappa_y = 10.^-2;

kappa_tau = 0.1; 
kappa_y = 0.5;


[wx, wy, wt, Ax_est, Ay_est, out] = mspoc(X_tr, Y_tr, ...
    'tau_vector', 0:(length(hrf)-1), ...
    'kappa_tau', kappa_tau, 'kappa_y', kappa_y);

%% apply to data

sy_est = wy' * Y;
px_est = zeros(1,Ne);
for e=1:Ne
    px_est(e) = var(X(:,:,e) * wx);
end
wt_tmp = wt/sum(wt);
px_flt_est = filter(wt_tmp, 1, px_est);


corr_tr = corrcoef(sy_est(tr_idx), px_flt_est(tr_idx));
corr_tr = (corr_tr(1,2));

corr_te = corrcoef(sy_est(te_idx), px_flt_est(te_idx));
corr_te = (corr_te(1,2));

corr_pat_y = corrcoef(Ay(:,1), Ay_est);
corr_pat_y = (corr_pat_y(1,2));

corr_pat_x = corrcoef(Ax(:,1), Ax_est);
corr_pat_x = (corr_pat_x(1,2));

fprintf('corr_tr = %g, corr_te = %g\n', abs(corr_tr), abs(corr_te))
fprintf('corr x pattern = %g , corr y pattern = %g \n', abs(corr_pat_x), abs(corr_pat_y));

%% plot results
figure,
rows = 4;
cols = 1;

subplot(rows,cols,1)
plot([zscore(px_flt_est)', sign(corr_tr)*zscore(sy_est)'])
title('estimated source time courses')

if Nt > 0
    subplot(rows,cols,2)
    sgn = sign(hrf*wt);
    % plot(zscore([hrf', sgn*wt, sgn*out.Atau]));
    plot(zscore([hrf', sgn*wt]));
    title('hrf vs wt')
end

subplot(rows,cols,3)
sgn = sign(Ax(:,1)'*Ax_est);
plot(zscore([Ax(:,1), sgn*Ax_est]))
title('Ax vs estimated Ax')

subplot(rows,cols,4)
sgn = sign(Ay(:,1)'*Ay_est);
plot(zscore([Ay(:,1), sgn*Ay_est]))
title('Ay vs estimated Ay')


