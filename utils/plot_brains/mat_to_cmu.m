function[data, meta] = mat_to_cmu(mat)
%MAT_TO_CMU  convert image(s) in matrix format to CMU format
%
% Usage: [data, meta] = mat_to_cmu(mat)
%
% INPUTS:
%          mat: a 4D matrix of brain images.
%
% OUTPUTS:
%         data: a cell array containing vectors of voxel activations.
%
%         meta: a struct with the following fields:
%            nvoxels: total number of voxels containing brain
%         coordToCol: dimx by dimy by dimz matrix of voxel numbers (zeros
%                     indicate no voxel at the corresponding location)
%         colToCoord: nvoxels by 3 matrix of voxel locations
%
%          
%
% SEE ALSO: CMU_TO_MAT, CONSTRUCT_META, PLOT_BRAIN2D, PLOT_BRAIN3D, SLICES, CAT
%
%  AUTHOR: Jeremy R. Manning
% CONTACT: manning3@princeton.edu

% CHANGELOG:
% 12-11-13 jrm  wrote it.

sliced_mats = slices(mat, ndims(mat));
meta = construct_meta(size(sliced_mats{1}));
data = cellfun(@(x)(x(:)'), sliced_mats, 'UniformOutput', false);

good_inds = ~isnan(mean(mat, ndims(mat)));
good_inds = good_inds(:)';

meta = meta_select_voxels(meta, find(good_inds));
data = cellfun(@(x)(x(good_inds)), data, 'UniformOutput', false);

