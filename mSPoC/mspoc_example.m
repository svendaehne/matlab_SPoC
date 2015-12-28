

%% params


n_channels = [30,50]; % nr of channels for X and for Y data
K = 1; % number of "target source pairs", i.e. sources with coupled envelopes
% n_sources = [30,400];
n_sources = n_channels;

Nt = 10;

band = [8,12]; % frequency band of interest here
SNR = 0.5; % signal-to-noise ratio (between 0 and 1) in terms of variance explained by the target source

Ne_tr = 500; % number of training epochs
Ne_te = 500; % number of test epochs
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
    hrf = hrf/sum(hrf);
    
    figure
    plot(hrf)
    xlabel('time lags')
    title('true transfer function (HRF)')
end


fprintf('\n')
fprintf('--- mSPoC example ---\n')
fprintf('\n')
fprintf('Number of channels in X dataset = %d\n', n_channels(1))
fprintf('Number of channels in Y dataset = %d\n', n_channels(2))
fprintf('Number of sources in X dataset = %d\n', n_sources(1))
fprintf('Number of sources in Y dataset = %d\n', n_sources(2))
fprintf('Number of shared sources = %d\n', K)



%% create example toy data

fprintf('Simulating data ... ')
[X, Y, Sx, Sx_env, Sy, Ax, Ay] = create_mspoc_example_data(K, n_sources, ...
            n_channels, hrf, SNR, samples_per_second, Ne, Te, band);

X = permute(reshape(X, [Te, Ne, n_channels(1)]), [1,3,2]);   
Y = Y';

Y_tr = Y(:,tr_idx);
X_tr = X(:,:,tr_idx);
fprintf(' done\n ')

%% set mspoc options

mspoc_params = [];
mspoc_params.tau_vector = 0:(length(hrf)-1);
mspoc_params.pca_Y_var_expl = 0.99;


%% optimize regularizers on training data

kappa_tau_list = 10.^(-1:1);
kappa_y_list = 10.^(-2:2);

% kappa_tau_list = 10.^0;
% kappa_y_list = 10.^0;

[best_kappa_tau, best_kappa_y, xval_out] = ...
    optimize_mspoc_regularizers(X_tr, Y_tr, mspoc_params ...
        ,'kappa_tau_list', kappa_tau_list ...
        ,'kappa_y_list', kappa_y_list ...
    );

%% mkSPoC on training data

mspoc_params.kappa_tau = best_kappa_tau; 
mspoc_params.kappa_y = best_kappa_y;
mspoc_params.verbose = 2;

[wx, wy, wt, Ax_est, Ay_est, out] = mspoc(X_tr, Y_tr, mspoc_params);

%% apply to training and test data

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
rows = 2;
cols = 3;

subplot(rows,cols,1:cols)
hold on
plot([zscore(px_flt_est)', sign(corr_tr)*zscore(sy_est)'])
title('estimated source time courses')
legend({'convolved power of estimated X source', 'time-course of estimated Y source'})
box on
plot(Ne_tr*[1,1], ylim, 'color',0.6*[1,1,1])
text(Ne_tr*1.01, 3, '--> test data')
h=text(Ne_tr*0.99, 3, 'training data <--');
set(h, 'HorizontalAlignment', 'right');


subplot(rows,cols,cols+1)
sgn = sign(Ax(:,1)'*Ax_est);
plot(zscore([Ax(:,1), sgn*Ax_est]))
title('spatial pattern of X source')
legend({'true a_x', 'estimated a_x'})


if Nt > 0
    subplot(rows,cols,cols+2)
    sgn = sign(hrf*wt);
%     plot(zscore([hrf', sgn*wt, sgn*out.Atau]));
    plot(zscore([hrf', sgn*wt]));
    title('hrf vs wt')
    legend({'true hrf', 'w_{\tau}'})
end


subplot(rows,cols,cols+3)
sgn = sign(Ay(:,1)'*Ay_est);
plot(zscore([Ay(:,1), sgn*Ay_est]))
title('spatial pattern of Y source')
legend({'true a_y', 'estimated a_y'})


