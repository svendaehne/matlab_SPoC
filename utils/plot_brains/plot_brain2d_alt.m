function[hs] = plot_brain2d_alt(x,meta,varargin)
%PLOT_BRAIN2D  Plot slices of a brain image
%
% Usage: plot_brain2d_alt(x,meta,[nrow],[ncol],[dim],[clim],[in_roi])
%
% INPUTS:
%      x: a 1 by nvoxels vector of voxel activations
%
%   meta: a struct with the following fields:
%        nvoxels: total number of voxels containing brain
%     coordToCol: dimx by dimy by dimz matrix of voxel numbers (zeros
%                 indicate no voxel at the corresponding location)
%     colToCoord: nvoxels by 3 matrix of voxel locations
%
%     **TIP: meta can also be an nvoxels by 3 matrix of voxel locations**
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
% OUTPUTS: [none]
%
% SEE ALSO: PLOT_BRAIN2D_ALT, PLOT_BRAIN3D, SANEPCOLOR, SLICES,
%           GETTIGHTSUBPLOTHANDLES, IMCONTOUR, SUBPLOT, PCOLOR,
%           IMAGESC, JET
%
%  AUTHOR: Jeremy R. Manning
% CONTACT: manning3@princeton.edu

% CHANGELOG:
% 4-11-12  jrm  wrote it.
% 2-22-13  jrm  minor edits to comments.
% 11-2-13  jrm  re-write as wrapper to plot_brain2d
% 12-12-13 jrm  fixed some typos in documentation
% 12-12-13 jrm  removed in_roi option.

img = cmu_to_mat(x, meta);
if isempty(varargin)
    hs = plot_brain2d(img);
else
    hs = plot_brain2d(img, varargin{:});
end