function[] = plot_brain3d(x,erode_layers)
%PLOT_BRAIN3D  Plot a 3d brain image
%
% Each voxel is plotted as a colored cube, such that the entire image forms
% a 3d brain.
%
% Usage: plot_brain3d(x,erode_layers)
%
% INPUTS:
%            x: a 3d matrix (tensor) of voxel activations
%      
% erode_layers: optional argument specifying how many layers to make much
%               more transparent.  this allows patterns beneath the surface
%               of the brain to be more-easily seen.
%
% OUTPUTS: [none]
%
% SEE ALSO: PLOT_BRAIN2D, PLOT_BRAIN3D, PATCH_3DARRAY, IMERODE
%
%  AUTHOR: Jeremy R. Manning
% CONTACT: manning3@princeton.edu

% CHANGELOG:
% 4-11-12  jrm  wrote it.
% 2-22-13  jrm  ensure voxels appear as cubes.
% 12-11-13 jrm  renamed to plot_brain3d

if ~exist('erode_layers','var')
    erode_layers = 0;
end

se = strel('arbitrary',ones([3 3 3]));
mask = ~isnan(x);
for i = 1:erode_layers
    mask = imerode(mask,se);
end

outer = x;
outer(mask) = nan;
inner = x;
inner(~mask) = nan;

hold on;
h1 = PATCH_3Darray(inner,'col');
if iscell(h1)
    for i = 1:length(h1)
        set(h1{i},'FaceAlpha',0.25);%,'EdgeColor','none');
    end
else
    set(h1,'FaceAlpha',0.25);
end

if erode_layers > 0
    h2 = PATCH_3Darray(outer,'col');
    if iscell(h2)
        for i = 1:length(h2)
            set(h2{i},'FaceAlpha',0.1,'EdgeColor','none');
        end
    else
        set(h2,'FaceAlpha',0.1,'EdgeColor','none');
    end
end

axis square;
axis off;

set(gca,'CameraPosition',[-239.1516 -192.6841 256.3707]);
