function [Y_w, My, Cyy] = util_mspoc__prepare_Y_signal(Y, varargin)

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, ...
    'pca_Y_var_expl', 0.95, ...
    'Cyy', []);


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
