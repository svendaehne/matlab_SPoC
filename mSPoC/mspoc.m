function [Wx, Wy, Wtau, Ax, Ay, out] = mspoc(X, Y, varargin)
% multimodal (temporal) (kernel) Source Power Co-modulation Analysis (mSPoC)
%
% Finds spatial filters wx and wy and a temporal filter wt such that the
% filtered power of the projected x-signal maximally covaries with the
% projected y-signal.
%
% [Wx, Wy, Wtau, Ax, Ay, lambda corr_values, p_values] = mspoc(X, Y, varargin)
%
% Input:
% size(X) is [n_samples_per_epoch_x, n_channels_x, n_epochs]
% size(Y) is [n_channels_y, n_epochs]
%   - n_samples_per_epoch_x should be enough to estimate power, i.e. at least
%   a few oscillations
%
% For more options see the first few lines of the code....
%
% Output:
% Wx - set of x-signal filters (columns)
% Wy - set of y-signal filters (columns)
% Wtau - set time filters (columns)
% Ax - set of x-signal patterns (in the columns)
% Ay - set of y-signal patterns (in the columns)
% lambda - vector of objective function values, i.e. covariance between the
%  components in X and Y
% corr_values - vector of correlations between the X and Y components
% p_values - p-values corresponding to corr_values
%
% sven.daehne@tu-berlin.de, 10-2011

%% params
opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, ...
    'n_component_sets', 1, ... % the number of components to extract
    'tau_vector', 0, ... % maximum timeshift of X relative to Y, given in # epochs
    'use_log', 0, ...
    'n_random_initializations', 10, ...
    'max_optimization_iterations', 20, ...
    'kappa_tau', 10^-3, ...
    'kappa_y', 10^-3, ...
    'pca_Y_var_expl', 0.95, ...
    'eps', 10^-3, ... % convergence criterium, applied to correlation in each iteration
    'mask', [], ...
    'Cxx', [], ...
    'Cxxe', [], ... trialwise covariance matrices of the X-signal
    'Cyy', [], ...
    'verbose',1, ...
    'component_idx', 1, ... % used internally, just for output
    'is_white_y', 0 ...
    );

if opt.verbose > 0
    fprintf('\n\n-> starting mSPoC analysis\n')
end

tau = opt.tau_vector;
Nt = length(tau);
if not(isempty(X))
    Nex = size(X,3);
elseif not(isempty(opt.Cxxe))
    Nex = size(opt.Cxxe,3);
else
    error('X or Cxxe must be given!')
end
Ney = size(Y,2);

if not(Nex == Ney)
    error('X and Y must have the same number of epochs!!!');
end
Ne = Nex;

if isempty(opt.mask)
    opt.mask = 0 < ones(1,Ne);
end
opt.mask(1:(length(tau)-1)) = 0;


%% whiten the X and Y signal and store the whitening matrices
if opt.verbose > 1
    fprintf('   start whitening\n')
end

[Cxxe_w, Cxxe, Cxx, Mx] = prepare_X_signal(X, opt);

% subtract temporal mean
Y = bsxfun(@minus, Y, mean(Y,2));
if isfield(opt, 'My')
    My = opt.My;
    Y_w = My'*Y;
else
    [Y_w, My] = prepare_Y_signal(Y, opt);
end
    

if opt.verbose > 1
    fprintf('   using %d X-components\n',size(Mx,2));
    fprintf('   using %d Y-components\n',size(My,2));
end

%% number of dims in X and Y, and number of samples per epoch in X
n_component_sets = opt.n_component_sets;
min_n = min(size(Mx,2), size(My,2));
if n_component_sets > min_n;
    n_component_sets = min_n;
end

%% remove outlier trials
x_rm_idx = detect_max_eigenvalue_epochs(Cxxe_w,[],0,2);

thr = mean(mean(Y_w))+2.5*std(mean(Y_w));
y_rm_idx = detect_high_variance_epochs(Y_w, thr,0);

rm_idx = union(x_rm_idx, y_rm_idx);
if opt.verbose > 1
    fprintf('   removing %d outlier trials.\n',length(rm_idx));
end
Cxxe_w(:,:,rm_idx) = [];
Y_w(:,rm_idx) = [];
opt.mask(rm_idx) = [];


%% optimize the filter sets (iteratively if more than one are requested)
Wx = zeros(size(Mx,2), n_component_sets);
Wy = zeros(size(My,2), n_component_sets);
Wtau = zeros(Nt, n_component_sets);
r_values = zeros(1, n_component_sets);

for k=1:n_component_sets
    if opt.verbose > 0
        fprintf('  optimizing component set %d/%d\n', k, n_component_sets)
    end
    
    % prepare deflation matrices for X and Y
    if k>1 % if we already have filers...
        % ...create a basis of the space that is orthogonal to them
        Bx = null(Wx');
        By = null(Wy');
    else
        % if k==1, no further projection is needed
        Bx = eye(size(Wx,1));
        By = eye(size(Wy,1));
    end
    
    % project whitended data into "deflated space"
    Cxxe_dfl = zeros(size(Bx,2), size(Bx,2), size(Cxxe_w,3));
    for n=1:size(Cxxe_w,3)
        Cxxe_dfl(:,:,n) = Bx' * Cxxe_w(:,:,n) * Bx;
    end
    Y_dfl = By' * Y_w;
    
    % perform mSPoC optimization in whitened and deflated space
    [wx, wy, wtau, r] = optimize_filters(Cxxe_dfl, Y_dfl, Mx*Bx, My*By, opt);
    
    % project weight vectors back into undeflated space
    Wx(:,k) = Bx * wx;
    Wy(:,k) = By * wy;
    Wtau(:,k) = wtau;
    r_values(k) = r;
end


%% create patterns

% project filers and patterns back into original sensor space
Wx = Mx * Wx;
Wy = My * Wy;

for k=1:size(Wx,2)
    Wx(:,k) = Wx(:,k) ./ sqrt(Wx(:,k)' * Cxx * Wx(:,k));
    Wy(:,k) = Wy(:,k) ./ std(Wy(:,k)' * Y);
end

Ax = Cxx * Wx / (Wx' * Cxx * Wx);
Sy = Wy' * Y;
Ay = Y * (Y' * Wy) / (Sy*Sy');

Nx = size(Cxxe,1);
Cxxe_vec = reshape(Cxxe, [Nx*Nx, Ne]);
corr_values = zeros(1,size(Ax,2));
Atau = zeros(size(Wtau));
for k=1:size(Ay,2)
    Wtau(:,k) = sign(sum(Wtau(:,k))) * Wtau(:,k);
    
    [~, mm_idx] = max(abs(Ax(:,k)));
    sgn = sign(Ax(mm_idx,k));
    Ax(:,k) = sgn * Ax(:,k);
    Wx(:,k) = sgn * Wx(:,k);
        
    sy = Wy(:,k)' * Y;
    [~, mm_idx] = max(abs(Ay(:,k)));
    c = sign(Y(mm_idx,:) * sy');
    Ay(:,k) = c * Ay(:,k);
    Wy(:,k) = c * Wy(:,k);
    
    wx_vec = reshape(Wx(:,k)*Wx(:,k)', [Nx*Nx,1]);
    px = wx_vec'*Cxxe_vec;
    px = px-mean(px);
    
    pxf = filter(Wtau(:,k), 1, px);
    tmp = corrcoef(pxf', sy');
    corr_values(k) = tmp(1,2);
    
    
    Pxe = zeros(Nt, Ne);
    for kk=1:Nt
        Pxe(kk,:) = circshift(px, [1, tau(kk)]);
    end
    Cpp = Pxe*Pxe';
    Atau(:,k) = Cpp * Wtau(:,k);
end


out = [];
out.corr_values = corr_values;
out.Atau = Atau;

% out.p_values = p_values;
% 
% out.W_spoc = Mx * out_aux.W_spoc;
% out.A_spoc = Cxx * out.W_spoc;
% out.ev_spoc = out_aux.ev_spoc;


function [wx, wy, wt, max_corr, aux_tmp] = optimize_filters(Cxxe, Y, BMx, BMy, opt)
% compute a single component pair in X and Y

tau = opt.tau_vector;
Nx = size(Cxxe,1);
Ny = size(Y,1);
Ne = size(Cxxe,3);
Nt = length(tau);

r_tmp = zeros(1,opt.n_random_initializations);
Wx_tmp = zeros(Nx, opt.n_random_initializations);
Wy_tmp = zeros(Ny, opt.n_random_initializations);
Wt_tmp = zeros(Nt, opt.n_random_initializations);
% out_aux = cell(1,opt.n_random_initializations);

Cxxe_vec = reshape(Cxxe, [Nx*Nx, Ne]);
mask = opt.mask > 0;
Cxx = squeeze(mean(Cxxe(:,:,mask),3));

kappa_y = opt.kappa_y;
kappa_tau = opt.kappa_tau;

Dyy = BMy'*BMy;
Dyy = Dyy / trace(Dyy);
Dpp = eye(Nt);
Cyy = cov(Y(:,mask)'); 
Cyy = Cyy / trace(Cyy);
Cyy_inv = inv(Cyy + kappa_y*Dyy);

% use an iterative approach, i.e. init wx randomly and
% iterate SPoC and CCA
for n=1:opt.n_random_initializations
    
    % power of projected x-signal
    wx = randn(Nx, 1);
    wx_vec = reshape(wx*wx', [Nx*Nx,1]);
    px = wx_vec'*Cxxe_vec;
    if opt.use_log
        px = log(px);
    end
        
    %% optimization loop
    ii = 0;
    converged = false;
    r_last_iter = -inf;
    r_iter = zeros(1,opt.max_optimization_iterations);
    Pxe = zeros(Nt, Ne);
    while ii < opt.max_optimization_iterations && not(converged)
        ii = ii + 1;
        
        %% use current wx to get ingredients for wt and beta_y
        % temporally embedded (zero-mean-) power signal
        px = px-mean(px(mask));
        for k=1:Nt
            Pxe(k,:) = circshift(px, [1, tau(k)]);
        end
        % CCA to get wt and beta_y
%         [wt, wy] = my_reg_CCA(Pxe(:,mask), Y(:,mask), kappa_tau, kappa_y, Dpp, Dyy, Cyy);
        [wt, wy] = my_reg_CCA2(Pxe(:,mask), Y(:,mask), kappa_tau, Dpp, Cyy_inv);
        
        
        %% use current wt and wy to get ingredients for new wx
        % filtered cov-time series
        Cxxe_vec_flt = filter(wt, 1, Cxxe_vec(:,mask)')';
        % projected y-signal, now univariate
        sy = wy' * Y;
        sy = (sy - mean(sy(mask)))/std(sy(mask));
        
        % apply SPoC to get wx
        Cxxz = reshape(sy(:,mask)*Cxxe_vec_flt', [Nx, Nx]) / Ne;
        [W, D] = eig(Cxxz, Cxx);
        [~, sorted_idx] = sort(abs(diag(D)), 'descend');
        W = W(:, sorted_idx);
        wx = W(:,1);
        
        % recompute (log-)power of projected x-signal
        wx_vec = reshape(wx*wx', [Nx*Nx,1]);
        px = wx_vec'*Cxxe_vec;
        if opt.use_log
            px = log(px);
        end
        px_flt = filter(wt, 1, px(mask)')';
        
        %% current correlation
        tmp_r = corrcoef([sy(mask)', px_flt']);
        r = tmp_r(1,2);
        r_iter(ii) = r;
        if opt.verbose > 2
            fprintf('   iter %02d --> abs(r) = %.3f \n', ii, abs(r))
        end
        if abs( abs(r)-abs(r_last_iter)) < opt.eps
            converged = true;
        end
        r_last_iter = r;
    end
%     r_iter(ii+1:end) = [];
    
    %% store correlation, wx, wy, and wt obtained in this run
    Wx_tmp(:,n) = wx;
    Wy_tmp(:,n) = wy;
    Wt_tmp(:,n) = wt;
    r_tmp(n) = r_last_iter;
    
    if opt.verbose
        fprintf('   Run %d done with corr = %.3f -> ', n, r_tmp(n))
        if not(converged)
            fprintf(' Not converged after %d iterations!\n', ii);
        else
            fprintf(' Converged after %d iterations!\n', ii);
        end
    end
end
% collect the best filter triple
[max_corr, max_idx] = max(abs(r_tmp));
wx = Wx_tmp(:,max_idx);
wy = Wy_tmp(:,max_idx);
wt = Wt_tmp(:,max_idx);
aux_tmp = [];

if opt.verbose
    fprintf('  Choosing result from run %d, with corr = %g\n', max_idx, max_corr);
end




function [wt, wy,r] = my_reg_CCA(Pxe, Ye, kappa_tau, kappa_y, Dp, Dy, Cyy)
[Nt,Ne] = size(Pxe);
Ny = size(Ye,1);
if isempty(Dp)
    Dp = eye(Nt);
end
if isempty(Dy)
    Dy = eye(Ny);
end

Pxe = bsxfun(@minus, Pxe, mean(Pxe,2));
Pxe = Pxe ./ (mean(std(Pxe,0,2)));

Cpp = Pxe * Pxe' / Ne;
% Cyy = Ye * Ye' / Ne;
Cpy = Pxe * Ye' / Ne;
Cyp = Ye * Pxe' / Ne;

% regularize the auto-covariance matrices
Cpp = Cpp/trace(Cpp) + kappa_tau*Dp/trace(Dp);
Cyy = Cyy/trace(Cyy) + kappa_y*Dy/trace(Dy);

A = [zeros(Nt), Cpy; Cyp, zeros(Ny)];
B = [Cpp, zeros(Nt,Ny); zeros(Ny,Nt), Cyy];

eig_opt.disp = 0;
[V,r] = eigs(A,B, 1, 'LM',eig_opt);
wt = V(1:Nt,1);
wy = V((Nt+1):end, 1);


function [wt, wy] = my_reg_CCA2(Pxe, Ye, kappa_tau, Dp, Cyy_inv)

Ne = size(Ye,2);

Pxe = bsxfun(@minus, Pxe, mean(Pxe,2));
Pxe = Pxe ./ (mean(std(Pxe,0,2)));

Cpp = Pxe * Pxe' / Ne;
% Cyy = Ye * Ye' / Ne;
Cpy = Pxe * Ye' / Ne;
Cyp = Ye * Pxe' / Ne;

% regularize the auto-covariance matrices
Cpp = Cpp/trace(Cpp) + kappa_tau*Dp/trace(Dp);

eig_opt.disp = 0;
Cyy_inv_Cyp = Cyy_inv * Cyp;
if size(Pxe,1) > 1
    M = (Cpp \ Cpy) * Cyy_inv_Cyp;
    [wt,r] = eigs(M, 1, 'LM',eig_opt);
else
    wt = 1;
end
wy = Cyy_inv_Cyp * wt;



function [Cxxe_w, Cxxe, Cxx, Mx] = prepare_X_signal(X, opt)
% spatially whiten the X signal and store the whitening matrix
if not(isempty(opt.Cxxe))
    Cxxe = opt.Cxxe;
else
    Nx = size(X,2);
    Ne = size(X,3);
    Cxxe = zeros(Nx, Nx, Ne);
    for e=1:Ne
        Cxxe(:,:,e) = cov(X(:,:,e));
    end
end

if not(isempty(opt.Cxx))
    Cxx = opt.Cxx;
else
    Cxx = squeeze(mean(Cxxe(:,:,opt.mask),3));
end

[V,D] = eig(Cxx);
[ev_sorted, sort_idx] = sort(diag(D), 'descend');

% compute an estimate of the rank of the data
tol = ev_sorted(1) * 10^-5;
r = sum(ev_sorted > tol);
n_components = r;

V = V(:,sort_idx);
Mx = V * diag(ev_sorted.^-0.5); % whitening filters are in the columns
Mx = Mx(:,1:n_components);

Cxxe_w = zeros(n_components,n_components,size(Cxxe,3));
for e=1:size(Cxxe,3)
    Cxxe_w(:,:,e) = Mx' * Cxxe(:,:,e) * Mx;
end


function [Y_w, My] = prepare_Y_signal(Y, opt)

% % PCA via singular value decomposition 
% [U,S] = svd(Y,'econ');
% ev_sorted = diag(S).^2; % PCA eigenvalues
% var_expl = cumsum(ev_sorted)./sum(ev_sorted);
% min_var_expl = opt.pca_Y_var_expl; 
% n = find(var_expl >= min_var_expl, 1, 'first');
% 
% % whitening matrix
% My = U * diag(1./diag(S));
% My = My(:,1:n);
% % whitened signal
% Y_w = My' * Y;


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
