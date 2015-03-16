function [W, A, r_values, all_r_values] = cspoc(X, maxmin_flag, varargin)
% canonical Source Power Correlations analysis (cSPoC)
% cSPoC extracts components with maximal envelope correlations from N oscillatory
% and multivariate datasets such as EEG/MEG/LFP/ECoG. For N>2, the
% extracted components maximize the pairwise averaged envelope
% correlations.
% 
% [W, A, r_values, all_r_values] = cspoc(X, maxmin_flag, varargin)
%
%
% Input:
%
% X             - Cell array containing the N datasets. The datasets can have
%                   a different number of channels but they must have a common
%                   time structure, i.e. they must have the same number of
%                   samples or epochs. 
%                   This means either size(X{n}) = [T, n_channels(n)] or
%                   size(X{n}) = [T(n), n_channels(n), n_epos]. So if the
%                   data is not epoched, then all datasets must have the
%                   same number of samples T but the number of channels can
%                   vary between datasets. If the data comes in epochs (or
%                   trials or time windows), i.e. X is 3-dimensional, then
%                   all datasets must have the same number of epochs.
% maxmin_flag   - Maximize (1) or minimize (-1) the envelope correlations
%
% Optional input arguments:
%
% 'n_component_sets' - Number of envelope-correlated components per dataset
%                       to be extracted. If the number of datasets is N=2,
%                       then we speak of a component pair. 
%                       Default: 1, i.e. extract only the highest correlating 
%                       component set/pair
% 'use_log'           - Optimize correlations of log-envelopes rather then
%                       envelopes. 
%                       Default: false
% 'average_over_epochs' - When optimizing the correlations, average the
%                       source envelopes within epochs. 
%                       Default: false
% 'n_repeats' - number of re-starts per component pair. Default: 10
% 'maxIter'   - maximum number of optimizer iterations. Default: 200
% 'pca_var_explained' - Dimensionality reduction via PCA: variance explained by
%                          PCA (must be value between 0 and 1), 
%                       Default: 1 (i.e. no dim reduction)
%                       NOTE: It is recommended to not use this parameter
%                       but to do dimensionality reduction by means of
%                       Spatio-Spectral Decomposition (SSD)
%
% Optional input arguments are given as keyword/value pair, 
% e.g. [...] = cspoc(X, 1, 'n_component_sets', 3, 'use_log', 0, ...).
%
%
% Output:
%
% W - Cell array of weight vectors that extract the source time courses from the data.
%          size(W{n}) = [n_channels(n), n_component_sets], i.e. filers are
%          in the columns.
% A - Cell array of spatial patterns that correspond to the respective components. 
%           Useful for visualisation. 
%           size(A{n}) = [n_channels(n), n_component_sets], i.e. patterns
%           are in the columns.
% r_values - Correlation values obtained between the envelopes of 
%            the extracted components. If N>2, the we average pairwise envelope
%            correlations. size(r_values) = n_component_sets
% all_r_values - Pairwise envelope correlations for each component set.
%                   size(all_r_values) = [N, N, n_component_sets]
%
%
% References:
% 
% S. Dahne, V. V. Nikulin, D. Ramirez, P. J. Schreier, K. R. Muller, S. Haufe, 
% "Finding brain oscillations with power dependencies in neuroimaging data", NeuroImage, 96:334-348, 2014 
%
% S. Haufe, S. Dahne, V. V. Nikulin, 
% "Dimensionality reduction for the analysis of brain oscillations"
% NeuroImage 101:583-597 2014
%
% sven.daehne@tu-berlin.de, 2013
% stefan.haufe@tu-berlin.de, 2013

%% params
opt= propertylist2struct(varargin{:});
opt= set_defaults(opt, ...
    'n_component_sets', 1, ... % number of component pairs to be extracted
    'use_log', 0, ...
    'n_repeats', 10, ... % number of re-starts per component pair
    'average_over_epochs', 0, ...
    'verbose', 0, ...
    'pca_var_explained', 1, ... % dim reduction: var explained by PCA on X data (between 0 and 1)
    'maxIter', 200, ... % number of optimizer iterations
    'derivative_check', 'off', ... % for internal use
    'Tx', [],... % for internal use
    'w_init', []); % for internal use


%% preprocessing

N = length(X);
if isempty(opt.Tx);
    opt.Tx = zeros(1,N);
    for n=1:N
        opt.Tx(n) = size(X{n},1);
    end
end
for n=1:N
    opt.Mx(n) = size(X{n},2);
end


[C, X_white, M, n_components] = prepare_data(X, opt.pca_var_explained);

if opt.verbose > 0
    fprintf('number of dimensions = ')
    for k=1:N
        fprintf('%d ', n_components(k));
    end
    fprintf('\n')
end
if strcmp(opt.n_component_sets, 'all')
    opt.n_component_sets = min(n_components);
end

opt.n_components_grouped = n_components;


%% optimize the filter sets (iteratively if more than one are requested)
W = cell(1,N);
r_values = zeros(1,opt.n_component_sets);
all_r_values = zeros(N, N, opt.n_component_sets);
for k=1:opt.n_component_sets
    
    if opt.verbose > 0
        fprintf('  optimizing component set %d/%d\n', k, opt.n_component_sets)
    end
    
    B = cell(1,N);
    if k>1 % if we already have filers...
        for n=1:N
            % ...create a basis of the space that is orthogonal to them
            B{n} = null(W{n}');
        end
    else
        for n=1:N
            B{n} = eye(n_components(n));
        end
    end
    
    X_redux = cell(1,N);
    for n=1:N
        % project the data to the space orthogonal to previous weight vectors
        % and thereby reduce the dimensionality
        X_redux{n} = X_white{n} * B{n};
    end
    n_components_redux = n_components - k + 1;
    opt_redux = opt;
    
    [w, r] = optimize_w(X_redux, maxmin_flag, n_components_redux, opt_redux);
    if opt.verbose > 1
        fprintf('   --> best corr = %g\n', r);
    end
    
    % compute averaged pair-wise power correlations
    p = zeros(size(X_redux{1},1),N);
    for n=1:N
        p(:,n) = abs(X_redux{n}*w{n});
    end
    if opt.average_over_epochs
        N_epos = size(p,1)/opt.Tx(1); % here we assume all data sets have the same epo-length. This will change. 
        p = squeeze(mean(reshape(p, [opt.Tx(1), N_epos, N])));
    end
    if opt.use_log
        p = log(p);
    end
        
    all_r_values(:, :, k) = corrcoef(p);
    r_values(k) = (sum(reshape(all_r_values(:, :, k), [], 1)) - N) / (N^2-N);
    
    % project the filters back to the space of original dimensionality
    for n=1:N
        W{n} = [W{n} B{n} * w{n}];
    end
end


%% post processing


% project back to original (un-whitened) channel space
for n=1:N
    W{n} = M{n}*W{n}; 
end

% Normalize the filters such that the extracted components have unit
% variance.
for n=1:N
    for k=1:size(W{n},2)
        W{n}(:,k) = W{n}(:,k) / sqrt(squeeze(W{n}(:,k)'*C{n}*W{n}(:,k)));
    end
end

% patterns
A = cell(1,N);
for n=1:N
    A{n} = C{n} * W{n} / (W{n}'*C{n}*W{n});
end


%% --- helper functions ---

%% prepare_data
function [Cxx, X_white, M, n_components] = prepare_data(X, min_var_explained)

N = length(X);
Cxx = cell(1,N);
X_white = cell(1,N);
M = cell(1,N);
n_components = zeros(1,N);
for ii = 1:N
    % sanity checks and hilbert transform
    if not(isempty(X{ii}))
        if ndims(X{ii}) == 3
            if all(isreal(X{ii}))
                n_epos = size(X{ii},3);
                for k=1:n_epos
                    X{ii}(:,:,k) = hilbert(X{ii}(:,:,k));
                end
            end
            % if x has 3 dims, we assume those to be [time, channels, epochs]
            X{ii} = reshape(permute(X{ii},[1,3,2]), [size(X{ii},1)*size(X{ii},3), size(X{ii},2)]);
        else
            % if x has 2 dims, we assume those to be [time, channels]
            if all(isreal(X{ii}))
                X{ii} = hilbert(X{ii});
            end
        end
    else
        error('X must be non-empty!');
    end
    % now size(X) = [time, channels]!
    X{ii} = X{ii} - repmat(mean(X{ii}), size(X{ii},1), 1);
    Cxx{ii} = cov(real(X{ii}));
end

for ii = 1:N
    rk = rank(Cxx{ii});
    % whiten the data
    [V, D] = eig(Cxx{ii});
    [ev_sorted, sort_idx] = sort(diag(D), 'descend');
    V = V(:,sort_idx);
    D = diag(ev_sorted);
    var_explained = cumsum(ev_sorted)/sum(ev_sorted);
    nc =  find(var_explained>=min_var_explained, 1);
    nc = min(nc, rk);
    M{ii} = V * diag(diag(D).^-0.5); % whitening filters are in the columns!!!
    M{ii} = M{ii}(:, 1:nc);
    X_white{ii} = X{ii}*M{ii};
    n_components(ii) = nc;
end


%% optimize_wxwy
function [w, r] = optimize_w(X_a, maxmin_flag, n_components, opt)

w_tmp = zeros(sum(n_components), opt.n_repeats);
fval_tmp = zeros(1, opt.n_repeats);
N = length(X_a);
for n=1:N
    X_a{n} = X_a{n}.';
end
% repeat optimization for given number of restarts
for k=1:opt.n_repeats
    if opt.verbose > 1
        fprintf('    starting optimisation run %d/%d: ',k,opt.n_repeats)
    end
    
    % check if we have some given initial startpoint
    if k==1 && not(isempty(opt.w_init))
        w_start = opt.w_init(1:size(w_tmp,1));
    else
        w_start = randn(size(w_tmp,1), 1);
    end
    
    % set optimizer options
    minFunc_opt = struct('DerivativeCheck', opt.derivative_check, 'MaxIter', opt.maxIter, ...
        'Display', 'off', 'useComplex', 0, 'numDiff', 0, 'TolFun', 1e-02, 'TolX', 1e-05);
    
    % call the optimizer: minFunc
    [w_tmp(:,k), fval_tmp(k), ~, output] = ...
        minFunc(@fval_grad, w_start, minFunc_opt, X_a, maxmin_flag, ...
            n_components, opt.use_log, opt.average_over_epochs, opt.Tx);
    
    if opt.verbose > 1
        fprintf(' f_value=%g n_iter=%d\n',fval_tmp(k),output.iterations);
    end
end
% find the optimal solution among the repeats
[r, idx] = min(fval_tmp);
if maxmin_flag>0
    r = -r;
end
w = w_vec_to_cell(w_tmp(:,idx), cell(1,N), n_components);




%% fval_grad
function [fval, grad] = fval_grad(w_vec, X_a, maxmin_flag, n_components, use_log, average_over_epochs, Tx)
% Computes the objective function value and the
% gradient at the current w

N = length(X_a);
w_cell = w_vec_to_cell(w_vec, cell(1,N), n_components);

N_epos = size(X_a{1},2);
if average_over_epochs
    N_epos = size(X_a{1},2)/Tx(1);
end

wtx = cell(1,N);
phi_x = zeros(N_epos,N);
grad_phi_x_wrt_wx = cell(1,N);

for n=1:N
    N_x = n_components(n);
    
    wx = w_cell{n};
    wtx{n} = wx'*X_a{n};
    
    px = abs(wtx{n});
    
    %     wtx_tmp = repmat(wtx{n}, N_x, 1);
    %     grad_px_wrt_wx = (real(wtx_tmp) .* real(X_a{n}) + imag(wtx_tmp) .* imag(X_a{n})) ./ abs(wtx_tmp);
    
    grad_px_wrt_wx = real(bsxfun(@times, wtx{n}./px, conj(X_a{n})));
    
    
    if average_over_epochs
        px = mean(reshape(px, [Tx(n), N_epos]));
        grad_px_wrt_wx = squeeze(mean(reshape(grad_px_wrt_wx', [Tx(n), N_epos, N_x])))';
    end
    if use_log
        phi_x(:,n) = log(px);
        grad_phi_x_wrt_wx{n} = bsxfun(@times, 1./px, grad_px_wrt_wx);
    else
        phi_x(:,n) = px;
        grad_phi_x_wrt_wx{n} = grad_px_wrt_wx;
    end
end

T = size(phi_x,1);
phi_x = bsxfun(@minus, phi_x, mean(phi_x));

% matrix of pair-wise envelope covariances
C_phi = phi_x' * phi_x / T;
% objective function value -> average pair-wise correlations
r = 0;
for n=1:N
    for k=(n+1):N
        r = r + C_phi(n,k)/sqrt(C_phi(n,n) * C_phi(k,k));
    end
end
norm_constant = (N^2-N)/2;
fval = -1*sign(maxmin_flag) * r / norm_constant;

% gradient
grad_x = cell(1,N);
for n=1:N
    grad_x{n} = zeros(n_components(n),1);
    for k=1:N
        delta_phi_nk = phi_x(:,k)' - phi_x(:,n)' * C_phi(n,k)/C_phi(n,n);
        grad_tmp = grad_phi_x_wrt_wx{n} * delta_phi_nk' / (T*sqrt(C_phi(n,n)*C_phi(k,k)));
        grad_x{n} = grad_x{n} + grad_tmp;
    end
end
tmp = zeros(sum(n_components),1);
grad = w_cell_to_vec(tmp, grad_x, n_components);
grad = -1*sign(maxmin_flag) * grad / norm_constant;


%% helper functions
function w_vec = w_cell_to_vec(w_vec, w_cell, n_dims)
stop = 0;
for k=1:length(n_dims)
    start = stop + 1;
    stop = start + n_dims(k) - 1;
    w_vec(start:stop) = w_cell{k};
end

function w_cell = w_vec_to_cell(w_vec, w_cell, n_dims)
stop = 0;
for k=1:length(n_dims)
    start = stop + 1;
    stop = start + n_dims(k) - 1;
    w_cell{k} = w_vec(start:stop);
end


