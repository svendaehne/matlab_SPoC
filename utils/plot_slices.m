function plot_slices(M, rows, cols, slice_dim)

tic

n_slices = size(M, slice_dim);
n_plots = rows*cols;
slices = round(linspace(1,n_slices, n_plots));

M = permute(M, [slice_dim, setdiff([1,2,3],slice_dim)]);

% figure
clim = max(abs(M(:)))*[-1,1];
for k=1:n_plots
    subplot(rows,cols,k)
    m = squeeze(M(slices(k),:,:));
    imagesc(m, clim)
    title(sprintf('slice nr %d', slices(k)))
    set(gca, 'xtick',[])
    set(gca, 'ytick',[])
end