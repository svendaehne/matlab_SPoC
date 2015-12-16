function [Cxxe_w, Cxxe, Cxx, Mx] = util_mspoc__prepare_X_signal(X, varargin)

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, ...
    'Cxx', [], ...
    'Cxxe', [] ... trialwise covariance matrices of the X-signal
    );


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