function[meta] = construct_meta(dims)
%CONSTRUCT_META  create a CMU meta matrix with the given dimensions
%
% Usage: meta = construct_meta(dims)
%
% INPUTS:
%         dims: a vector of finite positive integers specifying the desired
%               length of each dimension.
%
% OUTPUTS:         
%         meta: a struct with the following fields:
%            nvoxels: total number of voxels containing brain
%         coordToCol: dimx by dimy by dimz matrix of voxel numbers (zeros
%                     indicate no voxel at the corresponding location)
%         colToCoord: nvoxels by 3 matrix of voxel locations          
%
% SEE ALSO: MAT_TO_CMU, CMU_TO_MAT, FULLFACT, SLICES
%
%  AUTHOR: Jeremy R. Manning
% CONTACT: manning3@princeton.edu

% CHANGELOG:
% 12-11-13 jrm  wrote it.

if size(dims,1) > 1
    colToCoord = dims;
    dims = cellfun(@(x)(length(unique(x))), slices(dims,2));    
else
    colToCoord = fullfact(dims);
end
meta = struct('nvoxels', prod(dims), 'coordToCol', zeros(dims), 'colToCoord', colToCoord);

locs = cellfun(@mat2str, slices(meta.colToCoord, 1), 'UniformOutput', false);
locs = cellfun(@(x)(strrep(x,'[','(')), locs, 'UniformOutput', false);
locs = cellfun(@(x)(strrep(x,']',')')), locs, 'UniformOutput', false);
locs = cellfun(@(x)(strrep(x,' ',',')), locs, 'UniformOutput', false);

cmd = 'meta.coordToCol%s = %d;';
cellfun(@(x,i)(eval(sprintf(cmd,x,i))), locs, slices(1:meta.nvoxels));