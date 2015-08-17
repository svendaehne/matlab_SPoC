function[hs] = plot_brain2d(img,nrow,ncol,dim,c)
%PLOT_BRAIN2D  Plot slices of a brain image (3d matrix)
%
% Usage: plot_brain2d_alt(x,[nrow],[ncol],[dim],[clim])
%
% INPUTS:
%      x: a 3d matrix (tensor) of voxel activations
%
%   nrow: an optional argument specifying the number of rows of brain
%         images (default: 4)
% 
%   ncol: an optional argument specifying the number of columns of brain
%         images (default: 4)
%
%    dim: optional argument specifying the plane along which the 3d brain
%         image will be sliced. default: 3 (coronal slices).  options: 1
%         (sagittal), 2 (horizontal), 3 (coronal).
%
%   clim: optional argument specifying upper and lower bounds of the color
%         axis.  minimum value is mapped to blue; upper is mapped to red.
%         user can also input a matrix or brain image here (e.g. to use as
%         a reference for the to-be-plotted image.)  if so, the lower 10%
%         and upper 90% activation values of the given image are used as
%         the color axis bounds.
%
% OUTPUTS:
%     hs: a vector of handles to the subplots
%
% SEE ALSO: PLOT_BRAIN2D, PLOT_BRAIN3D, SANEPCOLOR, SLICES,
%           GETTIGHTSUBPLOTHANDLES, IMCONTOUR, SUBPLOT, PCOLOR,
%           IMAGESC, JET
%
%  AUTHOR: Jeremy R. Manning
% CONTACT: manning3@princeton.edu

% CHANGELOG:
% 2-22-13  jrm  wrote it.
% 11-2-13  jrm  rename to plot_brain2d
% 12-12-13 jrm  removed in_roi option.

warning('off','MATLAB:contour:ConstantData');

if ~exist('nrow','var') || isempty(nrow)
    nrow = 4;
end
if ~exist('ncol','var') || isempty(ncol)
    ncol = 4;
end
if ~exist('dim','var') || isempty(dim)
    dim = 3;
end
clf;

select_img = img;
nslices = min(nrow*ncol,size(img,dim));
brain_slices = slices(img,dim);
select_slices = slices(select_img,dim);

slice_edges = linspace(1,size(img,dim),nslices+1);
which_slices = mean([slice_edges(1:end-1) ; slice_edges(2:end)],1);
which_slices = unique(round(which_slices));

hs = getTightSubplotHandles(0.01,0.01,0.01,nrow,ncol);

if ~exist('c','var') || isempty(c)
    c = [prctile(img(:),10) prctile(img(:),90)];
elseif ~((min(size(c)) == 1) && (max(size(c)) == 2)) %c is a set of images...
    c = prctile(flatten(c), [10 90]);
end
for i = 1:length(hs)
    axes(hs(i)); %#ok<LAXES>
    if i <= length(which_slices)
        plot_helper(flipud(squeeze(brain_slices{which_slices(i)})'),flipud(squeeze(select_slices{which_slices(i)})'),c);%,xdim,ydim);        
    else        
        axis off;
    end
end

function[] = plot_helper(img,select_img,c)
hold on;
sanePColor(img); 
try
    caxis(c); 
end
colormap jet;
se = strel('arbitrary',ones([1 1]));
[~,h] = imcontour(imdilate(~isnan(select_img),se),1);
set(h,'LineColor','k','LineWidth',2);
hold off;
axis tight;
axis off;

function[f] = flatten(x)
if iscell(x)
    f = zeros(1, sum(cellfun(@numel, x)));
    start = 1;
    for i = 1:length(x)
        f(start:(start + numel(x{i})-1)) = x{i}(:)';
        start = start + numel(x{i});
    end
else
    f = x(:)';
end
