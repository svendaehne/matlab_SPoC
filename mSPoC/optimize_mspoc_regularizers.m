function [best_kappa_tau, best_kappa_y, out] = ...
    optimize_mspoc_regularizers(X, Y, mspoc_params, varargin)

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt ...
    ,'n_xvalidation_folds', 5 ...
    ,'kappa_tau_list', 10.^(-4:1) ...
    ,'kappa_y_list', 10.^(-4:1) ...
    ,'Cxxe', [] ...
    );

kappa_tau_list = opt.kappa_tau_list;
kappa_y_list = opt.kappa_y_list;


if not( isempty(X))
    Nx = size(X,2);
    Ne = size(X,3);
    
    Cxxe = zeros(Nx, Nx, Ne);
    for e=1:Ne
        Cxxe(:,:,e) = cov(squeeze(X(:,:,e)));
    end
elseif not(isempty(opt.Cxxe))
    Cxxe = opt.Cxxe;
    Nx = size(Cxxe,1);
    Ne = size(Cxxe,3);
end

Ny = size(Y,1);
Nt = 1;
if isfield(mspoc_params, 'tau_vector')
    Nt = length(mspoc_params.tau_vector);
end

n_folds = opt.n_xvalidation_folds;
n_tr_epos = floor(Ne/n_folds);

Cxxe_vec = reshape(Cxxe, [Nx*Nx, Ne]);

n_kappa_tau = length(kappa_tau_list);
n_kappa_y = length(kappa_y_list);
corr_tr = zeros(n_kappa_tau, n_kappa_y, n_folds);
corr_te = zeros(n_kappa_tau, n_kappa_y, n_folds);

Ax_all = zeros(n_kappa_tau, n_kappa_y, n_folds, Nx);
Ay_all = zeros(n_kappa_tau, n_kappa_y, n_folds, Ny);
Wt_all = zeros(n_kappa_tau, n_kappa_y, n_folds, Nt);

pattern_similarity_x = zeros(n_kappa_tau, n_kappa_y);
pattern_similarity_y = zeros(n_kappa_tau, n_kappa_y);
pattern_similarity_tau = zeros(n_kappa_tau, n_kappa_y);

display('cross-validating mSPoC parameters')

%% precompute some useful stuff
fprintf(' precomputing whitening for Y   ')
tr_idx = cell(1,n_folds);
te_idx = cell(1,n_folds);
Cyy_tr = cell(1,n_folds);
My_tr = cell(1,n_folds);
for k=1:n_folds
    
    fprintf('.')
    % compute train and test indices
    start_idx = (k-1)*n_tr_epos + 1;
    stop_idx = k*n_tr_epos;
    te_idx{k} = start_idx:stop_idx;
    tr_idx{k} = setdiff(1:Ne, te_idx{k});
    
    Y_tr = Y(:,tr_idx{k});
    Y_tr = bsxfun(@minus, Y_tr, mean(Y_tr,2));
    
    mspoc_params.Cyy = [];
    [~, My_tr{k}, Cyy_tr{k}] = prepare_Y_signal(Y_tr, mspoc_params);
    
end
fprintf('  done\n')

%% xval loops
for ii=1:n_kappa_tau
    for jj=1:n_kappa_y
        fprintf(' starting xval for kappa_tau = %g, kappa_y = %g ',kappa_tau_list(ii),kappa_y_list(jj))
        for k=1:n_folds
            
            %% run mSPoC with the current set of regularizers
            mspoc_params.kappa_tau = kappa_tau_list(ii);
            mspoc_params.kappa_y = kappa_y_list(jj);
            mspoc_params.n_component_sets = 1;
            mspoc_params.verbose = 0;
            mspoc_params.Cxxe = Cxxe(:,:,tr_idx{k});
            mspoc_params.Cyy = Cyy_tr{k};
            mspoc_params.My = My_tr{k};
            
            [wx, wy, wt, Ax_est, Ay_est] = mspoc([], Y(:,tr_idx{k}), mspoc_params);
            
            Ax_all(ii,jj,k,:) = Ax_est(:,1);
            Ay_all(ii,jj,k,:) = Ay_est(:,1);
            Wt_all(ii,jj,k,:) = wt(:,1);
            
            
            %% apply weight vectors
            sy_est = wy' * Y;
            wx_vec = reshape(wx*wx',[Nx*Nx,1]);
            px_est = wx_vec' * Cxxe_vec;
            
            px_flt_est = filter(wt, 1, px_est);
            
            %% compute correlations on training and test data
            tmp = corrcoef(sy_est(tr_idx{k}), px_flt_est(tr_idx{k}));
            corr_tr(ii,jj,k) = abs(tmp(1,2));
            tmp = corrcoef(sy_est(te_idx{k}), px_flt_est(te_idx{k}));
            corr_te(ii,jj,k) = abs(tmp(1,2));
            
            fprintf('.')
        end
        fprintf('  done\n')
        fprintf(' -> corr_tr = %0.3f +/- %0.3f, corr_te = %0.3f +/- %0.3f\n', ...
            mean(corr_tr(ii,jj,:)), std(corr_tr(ii,jj,:)), ...
            mean(corr_te(ii,jj,:)), std(corr_te(ii,jj,:)))
        
        % compute pattern similarity
        pattern_similarity_x(ii,jj) = pattern_similarity(squeeze(Ax_all(ii,jj,:,:))');
        pattern_similarity_y(ii,jj) = pattern_similarity(squeeze(Ay_all(ii,jj,:,:))');
        pattern_similarity_tau(ii,jj) = pattern_similarity(squeeze(Wt_all(ii,jj,:,:))');
    end
end

avg_corr_te = mean(corr_te,3);
[max_corr, idx] = max(abs(avg_corr_te(:)));
[ii,jj] = ind2sub([n_kappa_tau, n_kappa_y], idx);

best_kappa_tau = kappa_tau_list(ii);
best_kappa_y = kappa_y_list(jj);

out = [];
out.corr_tr = corr_tr;
out.corr_te = corr_te;
out.kappa_tau_list = kappa_tau_list;
out.kappa_y_list = kappa_y_list;
out.max_corr = max_corr;

out.pattern_similarity_x = pattern_similarity_x;
out.pattern_similarity_y = pattern_similarity_y;
out.pattern_similarity_tau = pattern_similarity_tau;



function sim = pattern_similarity(P)
N = size(P,2);
S = P' * P;
d = sqrt(diag(S));
S = abs(S) ./ (d * d');
sim = (sum(S(:))-N) / (N*(N-1));

function [Y_w, My, Cyy] = prepare_Y_signal(Y, opt)

% compute covariance matrix
[Ny, Ne] = size(Y);
if isempty(opt.Cyy)
    if Ne > Ny
        % if there are more samples than dimensions, compute the spatial
        % covariance matrix
        Cyy = Y*Y' / Ne;
    else
        % if there are more dimensions than samples, compute the temporal
        % covariance matrix
        Cyy = Y'*Y;
    end
else
    Cyy = opt.Cyy;
end

[V,D] = eig(Cyy);
[ev_sorted, sort_idx] = sort(diag(D), 'descend');
V = V(:,sort_idx);
My = V * diag(ev_sorted.^-0.5); % whitening filters are in the columns
% PCA dim-reduction
var_expl = cumsum(ev_sorted)./sum(ev_sorted);
min_var_expl = opt.pca_Y_var_expl; 
n = find(var_expl >= min_var_expl, 1, 'first');
My = My(:,1:n);

% whitening and possible dim-reduction
if Ne > Ny
    Y_w = My' * Y;
else
    Y_w = diag(std(V(:,1:n))) \ V(:,1:n)';
    My = sqrt(Ne) * Y * (V(:,1:n) * diag(1./ev_sorted(1:n)'));
end

