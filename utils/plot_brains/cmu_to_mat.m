function[img] = cmu_to_mat(data, meta)
%CMU_TO_MAT  convert image(s) in CMU format to matrix format
%
% Usage: img = cmu_to_mat(data, meta);
%
% INPUTS:
%         data: either a matrix whose rows are vectors of voxel activations
%               or a cell array whose entries are row vectors of voxel
%               activations.
%
%         meta: a struct with the following fields:
%            nvoxels: total number of voxels containing brain
%         coordToCol: dimx by dimy by dimz matrix of voxel numbers (zeros
%                     indicate no voxel at the corresponding location)
%         colToCoord: nvoxels by 3 matrix of voxel locations
%
%     **TIP: meta can also be an nvoxels by 3 matrix of voxel locations**
%
% OUTPUTS:
%          img: a 4D matrix of brain images.  (if data is a vector, or if
%               length(data) = 1, then img will be a 3D matrix.
%
% SEE ALSO: MAT_TO_CMU, PLOT_BRAIN2D, PLOT_BRAIN3D, SLICES, CAT
%
%  AUTHOR: Jeremy R. Manning
% CONTACT: manning3@princeton.edu

% CHANGELOG:
% 12-11-13 jrm  wrote it.
% 12-12-13 jrm  also accept matrices of voxel activations and voxel
%               locations

if ~iscell(data)
    data = slices(data, 1);
end
if ~isstruct(meta)
    meta = struct('colToCoord', meta, 'nvoxels', length(data{1}), 'coordToCol', zeros(max(meta)));
end

img = cell(size(data));
template = nan(size(meta.coordToCol));

for i = 1:length(data)
    next = template;
    if ndims(template) == 2 %#ok<ISMAT>
        for j = 1:meta.nvoxels
            next(meta.colToCoord(j,1),meta.colToCoord(j,2)) = data{i}(j);
        end
    elseif ndims(template) == 3
        for j = 1:meta.nvoxels
            next(meta.colToCoord(j,1),meta.colToCoord(j,2),meta.colToCoord(j,3)) = data{i}(j);
        end
    else
        for j = 1:meta.nvoxels
            eval(sprintf('next(%s) = data{i}(j);',join(',',arrayfun(@num2str,meta.colToCoord(j,:),'UniformOutput',false))));
        end
    end    
    img{i} = next;
end
img = cat(ndims(template)+1, img{:});