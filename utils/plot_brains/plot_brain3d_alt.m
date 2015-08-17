function[] = plot_brain3d(x,meta,erode_layers)
%PLOT_BRAIN3D  Plot a 3d brain image
%
% Each voxel is plotted as a colored cube, such that the entire image forms
% a 3d brain.
%
% Usage: plot_brain3d(x,meta,erode_layers)
%
% INPUTS:
%            x: a 1 by nvoxels vector of voxel activations
%
%         meta: a struct with the following fields:
%            nvoxels: total number of voxels containing brain
%         coordToCol: dimx by dimy by dimz matrix of voxel numbers (zeros
%                     indicate no voxel at the corresponding location)
%         colToCoord: nvoxels by 3 matrix of voxel locations
%
%     **TIP: meta can also be an nvoxels by 3 matrix of voxel locations**
%      
% erode_layers: optional argument specifying how many layers to make much
%               more transparent.  this allows patterns beneath the surface
%               of the brain to be more-easily seen.
%
% OUTPUTS: [none]
%
% SEE ALSO: PLOT_BRAIN2D, PLOT_BRAIN3D_ALT, PATCH_3DARRAY, IMERODE
%
%  AUTHOR: Jeremy R. Manning
% CONTACT: manning3@princeton.edu

% CHANGELOG:
% 4-11-12  jrm  wrote it.
% 2-22-13  jrm  ensure voxels appear as cubes.
% 12-11-13 jrm  re-wrote as wrapper to plot_brain3d
% 12-12-13 jrm  clarified necessary fields in meta struct

if ~exist('erode_layers', 'var')
    erode_layers = 0;
end

img = cmu_to_mat(x, meta);
plot_brain3d(img, erode_layers);

