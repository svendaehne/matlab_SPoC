function [varargout] = PATCH_3Darray(varargin)
% PATCH_3Darray  Plot a 3D array using patch to create a quadrangular surface mesh
%==========================================================================
% FILENAME:          PATCH_3Darray.m
% AUTHOR:            Adam H. Aitkenhead
% INSTITUTION:       The Christie NHS Foundation Trust
% CONTACT:           adam.aitkenhead@physics.cr.man.ac.uk
% DATE:              17th August 2010
% PURPOSE:           Plot a 3D array using patch to create a quadrangular
%                    surface mesh
%
% USAGE:             [hpat] = PATCH_3Darray(gridINPUT)
%           ..or..   [hpat] = PATCH_3Darray(gridINPUT,cmap)
%           ..or..   [hpat] = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ)
%           ..or..   [hpat] = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,cmap)
%           ..or..   [hpat] = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,'col')
%           ..or..   [hpat] = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,'sym')
%           ..or..   [hpat] = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,cmap,'col')
%           ..or..   [hpat] = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,cmap,'sym')
%           ..or..   [hpat,hcbar] = PATCH_3Darray(gridINPUT,gridX,gridY,gridZ,cmap,'sym')
%
% INPUT PARAMETERS:
%   gridINPUT        - 3D array of size (P,Q,R)
%                      If not using the flags 'col' or 'sym', then
%                      gridINPUT should be a logical array.  If using the
%                      flags 'col' or 'sym', then gridINPUT should be a
%                      numeric array where voxels which are not to be
%                      displayed contain the value NaN. 
%   gridX (optional) - A 1xP array  - List of the X axis coordinates.
%   gridY (optional) - A 1xQ array  - List of the Y axis coordinates.
%   gridZ (optional) - A 1xR array  - List of the Z axis coordinates.
%   gridZ (optional) - A 1xR array  - List of the Z axis coordinates.
%   cmap  (optional) - A Nx3 array  - The colormap definition.  When
%                      plotting using the 'col' or 'sym' flag, cmap must be
%                      an Nx3 array, eg jet(32).  When plotting a logical
%                      array, cmap must be an RGB triplet, eg [0.5,0.5,0]. 
%
% ADDITIONAL INPUT FLAGS:
%   'col'  (optional) - When this flag is present, the colour of each
%                       facet corresponds to the value of the voxel.  In
%                       the input array gridINPUT, voxels which are not to 
%                       be displayed should have a value of NaN.
%   'sym'  (optional) - Similar to the flag 'col', except that the
%                       colorbar is symmetric about zero.
%   'barN' (optional) - Display a colorbar on North of plot.
%   'barE' (optional) - Display a colorbar on East of plot.
%   'barS' (optional) - Display a colorbar on South of plot.
%   'barW' (optional) - Display a colorbar on West of plot.
%
% OUTPUT PARAMETERS:
%   hpat  (optional)  - Handle to the patch object.
%   hcbar (optional)  - Handle to the colorbar.
%==========================================================================

%==========================================================================
% VERSION  USER  CHANGES
% -------  ----  -------
% 100817   AHA   Original version
% 101126   AHA   Better handling of objects which are only 1 pixel wide.
% 110310   AHA   Can now produce plots where the colour of each facet
%                corresponds to the value of the voxel.
% 110311   AHA   The user can now specify the colour map for the display of
%                a 3D numeric array.
% 110316   AHA   Set the facet edges to be a darker version of the face
%                colour.
% 110324   AHA   Bugfix:  Ensure the facet values lie within [1:cres]
%                rather than [0:cres].  Otherwise the 3D object will have a
%                hole, since any 0-value facets will not be plotted.
% 110401   AHA   Minor speed up by reducing number of calls to ind2sub.
% 110407   AHA   Added control of the colorbar position.
%==========================================================================


%======================================================
% CHECK THE INPUTS
%======================================================

checkchar = false(nargin,1);
checknl   = false(nargin,1);
for loopC = 1:nargin
  checkchar(loopC) = ischar(varargin{loopC});
  checknl(loopC)   = isnumeric(varargin{loopC}) + islogical(varargin{loopC});
end

%Check for numeric or logical inputs
if sum(checknl)==1
  checknlIND = find(checknl==1);
  gridINPUT  = varargin{checknlIND};
  gridSIZE   = size(gridINPUT);
  gridX      = [1:gridSIZE(1)];
  gridY      = [1:gridSIZE(2)];
  gridZ      = [1:gridSIZE(3)];
elseif sum(checknl)==2
  checknlIND = find(checknl==1);
  gridINPUT  = varargin{checknlIND(1)};
  gridSIZE   = size(gridINPUT);
  gridX      = [1:gridSIZE(1)];
  gridY      = [1:gridSIZE(2)];
  gridZ      = [1:gridSIZE(3)];
  cmap       = varargin{checknlIND(2)};
elseif sum(checknl)==4
  checknlIND = find(checknl==1);
  gridINPUT  = varargin{checknlIND(1)};
  gridSIZE   = size(gridINPUT);
  gridX      = varargin{checknlIND(2)};
  gridY      = varargin{checknlIND(3)};
  gridZ      = varargin{checknlIND(4)};
elseif sum(checknl)==5
  checknlIND = find(checknl==1);
  gridINPUT  = varargin{checknlIND(1)};
  gridSIZE   = size(gridINPUT);
  gridX      = varargin{checknlIND(2)};
  gridY      = varargin{checknlIND(3)};
  gridZ      = varargin{checknlIND(4)};
  cmap       = varargin{checknlIND(5)};
else
  error('Incorrect number of numeric/logical inputs.')
end

%Set defaults for colormap and colorbar options
symYN = 0;
colYN = 0;
cbar  = 0;
  
%Check for character inputs
if sum(checkchar)>=1
  checkcharIND = find(checkchar==1);
  for loopCH = 1:sum(checkchar)
    if strncmpi('sym',varargin{checkcharIND(loopCH)},3)==1
      symYN    = 1;
      colYN = 1;
    elseif strncmpi('col',varargin{checkcharIND(loopCH)},3)==1
      colYN = 1;
    elseif strcmpi('barno',varargin{checkcharIND(loopCH)})==1
      cbar = 0;
    elseif strcmpi('barN',varargin{checkcharIND(loopCH)})==1
      cbar = 1;
    elseif strcmpi('barE',varargin{checkcharIND(loopCH)})==1
      cbar = 2;
    elseif strcmpi('barS',varargin{checkcharIND(loopCH)})==1
      cbar = 3;
    elseif strcmpi('barW',varargin{checkcharIND(loopCH)})==1
      cbar = 4;
    end
  end
end

%Check that the colormap is a 1x3 array if a logical array is to be plotted
if colYN==0 && exist('cmap','var')==1 && numel(cmap)~=3
  error('When displaying a logical array, the colour map must be a 1x3 RGB triplet.')
end

%Check the size of the grid
if size(gridX,1)>size(gridX,2)
  gridX = gridX';
end
if size(gridY,1)>size(gridY,2)
  gridY = gridY';
end
if size(gridZ,1)>size(gridZ,2)
  gridZ = gridZ';
end

if ~isequal(gridSIZE,[numel(gridX),numel(gridY),numel(gridZ)])
  error('The dimensions of gridINPUT do not match the dimensions of gridX, gridY, gridZ.')
end


%======================================================
% OBTAIN A LOGICAL ARRAY FROM gridINPUT
%======================================================

if colYN == 1;
  gridINPUTlogical = ~isnan(gridINPUT);
elseif colYN == 0;
  gridINPUTlogical = gridINPUT;
  nanvoxels = isnan(gridINPUT);
  nancount  = sum(nanvoxels(:));
  if nancount>0
    gridINPUTlogical(~nanvoxels) = 1;
    gridINPUTlogical(nanvoxels)  = 0;
  end
  gridINPUTlogical = logical(gridINPUTlogical);
end


%======================================================
% REMOVE ANY OUTER UNUSED AREAS FROM gridINPUT
%======================================================

objectIND                 = find(gridINPUTlogical);
[objectX,objectY,objectZ] = ind2sub(gridSIZE,objectIND);

if objectX(1)~=objectX(end)
  gridINPUT        = gridINPUT(min(objectX):max(objectX),:,:);
  gridINPUTlogical = gridINPUTlogical(min(objectX):max(objectX),:,:);
  gridX            = gridX(min(objectX):max(objectX));
end
if objectY(1)~=objectY(end)
  gridINPUT        = gridINPUT(:,min(objectY):max(objectY),:);
  gridINPUTlogical = gridINPUTlogical(:,min(objectY):max(objectY),:);
  gridY            = gridY(min(objectY):max(objectY));
end
if objectZ(1)~=objectZ(end)
  gridINPUT        = gridINPUT(:,:,min(objectZ):max(objectZ));
  gridINPUTlogical = gridINPUTlogical(:,:,min(objectZ):max(objectZ));
  gridZ            = gridZ(min(objectZ):max(objectZ));
end

gridSIZE = size(gridINPUT);


%======================================================
% DEFINE THE LOWER AND UPPER LIMITS OF EACH VOXEL
%======================================================

if gridSIZE(1)==1
  gridXlower = gridX;
  gridXupper = gridX;
else
  gridXsteps = gridX(2:end)-gridX(1:end-1);
  gridXlower = gridX-[gridXsteps(1),gridXsteps]/2;
  gridXupper = gridX+[gridXsteps,gridXsteps(end)]/2;
end

if gridSIZE(2)==1
  gridYlower = gridY;
  gridYupper = gridY;
else
  gridYsteps = gridY(2:end)-gridY(1:end-1);
  gridYlower = gridY-[gridYsteps(1),gridYsteps]/2;
  gridYupper = gridY+[gridYsteps,gridYsteps(end)]/2;
end

if gridSIZE(3)==1
  gridZlower = gridZ;
  gridZupper = gridZ;
else
  gridZsteps = gridZ(2:end)-gridZ(1:end-1);
  gridZlower = gridZ-[gridZsteps(1),gridZsteps]/2;
  gridZupper = gridZ+[gridZsteps,gridZsteps(end)]/2;
end

%======================================================
% FOR EACH VOXEL, IDENTIFY WHETHER ITS 6 NEIGHBOURS ARE WITHIN THE OBJECT.
% IF ANY NEIGHBOUR IS OUTSIDE THE OBJECT, DRAW FACETS BETWEEN THE VOXEL AND
% THAT NEIGHBOUR.
%======================================================

gridINPUTshifted = false(gridSIZE);
if gridSIZE(1)>2
  gridINPUTwithborder = cat(1,false(1,gridSIZE(2),gridSIZE(3)),gridINPUTlogical,false(1,gridSIZE(2),gridSIZE(3)));         %Add border
  gridINPUTshifted    = cat(1,false(1,gridSIZE(2),gridSIZE(3)),gridINPUTshifted,false(1,gridSIZE(2),gridSIZE(3)));      %Add border
  gridINPUTshifted    = gridINPUTshifted + circshift(gridINPUTwithborder,[-1,0,0]) + circshift(gridINPUTwithborder,[1,0,0]);
  gridINPUTshifted    = gridINPUTshifted(2:end-1,:,:);  %Remove border
end
if gridSIZE(2)>2
  gridINPUTwithborder = cat(2,false(gridSIZE(1),1,gridSIZE(3)),gridINPUTlogical,false(gridSIZE(1),1,gridSIZE(3)));         %Add border
  gridINPUTshifted    = cat(2,false(gridSIZE(1),1,gridSIZE(3)),gridINPUTshifted,false(gridSIZE(1),1,gridSIZE(3)));      %Add border
  gridINPUTshifted    = gridINPUTshifted + circshift(gridINPUTwithborder,[0,-1,0]) + circshift(gridINPUTwithborder,[0,1,0]);
  gridINPUTshifted    = gridINPUTshifted(:,2:end-1,:);  %Remove border
end
if gridSIZE(3)>2
  gridINPUTwithborder = cat(3,false(gridSIZE(1),gridSIZE(2),1),gridINPUTlogical,false(gridSIZE(1),gridSIZE(2),1));         %Add border
  gridINPUTshifted    = cat(3,false(gridSIZE(1),gridSIZE(2),1),gridINPUTshifted,false(gridSIZE(1),gridSIZE(2),1));      %Add border
  gridINPUTshifted    = gridINPUTshifted + circshift(gridINPUTwithborder,[0,0,-1]) + circshift(gridINPUTwithborder,[0,0,1]);
  gridINPUTshifted    = gridINPUTshifted(:,:,2:end-1);  %Remove border
end

%Identify the voxels which are at the boundary of the object:
edgevoxellogical = gridINPUTlogical & gridINPUTshifted<6;
edgevoxelindices = find(edgevoxellogical)';
edgevoxelcount   = numel(edgevoxelindices);

%Calculate the number of facets there will be in the final quadrangular mesh:
facetcount = (edgevoxelcount*6 - sum(gridINPUTshifted(edgevoxelindices)) );

%Create an array to record...
%Cols 1-6: Whether each edge voxel's 6 neighbours are inside or outside the object.
neighbourlist = false(edgevoxelcount,6);

%Initialise arrays to store the quadrangular mesh data:
facetsALL  = zeros(facetcount,3,4);
normalsALL = zeros(facetcount,3);
valuesALL  = zeros(facetcount,1);

%Create a counter to keep track of how many facets have been written as the
%following 'for' loop progresses:
facetcountsofar = 0;

[subXALL,subYALL,subZALL] = ind2sub(gridSIZE,edgevoxelindices);

for loopP = 1:edgevoxelcount
  
  subX = subXALL(loopP);
  subY = subYALL(loopP);
  subZ = subZALL(loopP);
  
  if subX==1
    neighbourlist(loopP,1) = 0;
  else
    neighbourlist(loopP,1) = gridINPUTlogical(subX-1,subY,subZ);
  end
  if subY==1
    neighbourlist(loopP,2) = 0;
  else
    neighbourlist(loopP,2) = gridINPUTlogical(subX,subY-1,subZ);
  end
  if subZ==gridSIZE(3)
    neighbourlist(loopP,3) = 0;
  else
    neighbourlist(loopP,3) = gridINPUTlogical(subX,subY,subZ+1);
  end
  if subY==gridSIZE(2)
    neighbourlist(loopP,4) = 0;
  else
    neighbourlist(loopP,4) = gridINPUTlogical(subX,subY+1,subZ);
  end
  if subZ==1
    neighbourlist(loopP,5) = 0;
  else
    neighbourlist(loopP,5) = gridINPUTlogical(subX,subY,subZ-1);
  end
  if subX==gridSIZE(1)
    neighbourlist(loopP,6) = 0;
  else
    neighbourlist(loopP,6) = gridINPUTlogical(subX+1,subY,subZ);
  end
  
  facetCOtemp         = zeros(6-sum(neighbourlist(loopP,:)),3,4);
  normalCOtemp        = zeros(6-sum(neighbourlist(loopP,:)),3);
  facetcountthisvoxel = 0;
  
  if neighbourlist(loopP,1)==0   %Neighbouring voxel in the -x direction
    facetcountthisvoxel                    = facetcountthisvoxel+1;
    facetcountsofar                        = facetcountsofar+1;
    facetCOtemp(facetcountthisvoxel,1:3,1) = [ gridXlower(subX),gridYlower(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,2) = [ gridXlower(subX),gridYlower(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,3) = [ gridXlower(subX),gridYupper(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,4) = [ gridXlower(subX),gridYupper(subY),gridZlower(subZ) ];
    normalCOtemp(facetcountthisvoxel,1:3)  = [-1,0,0];
  end
  if neighbourlist(loopP,2)==0   %Neighbouring voxel in the -y direction
    facetcountthisvoxel                    = facetcountthisvoxel+1;
    facetcountsofar                        = facetcountsofar+1;
    facetCOtemp(facetcountthisvoxel,1:3,1) = [ gridXlower(subX),gridYlower(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,2) = [ gridXupper(subX),gridYlower(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,3) = [ gridXupper(subX),gridYlower(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,4) = [ gridXlower(subX),gridYlower(subY),gridZupper(subZ) ];
    normalCOtemp(facetcountthisvoxel,1:3)  = [0,-1,0];
  end
  if neighbourlist(loopP,3)==0   %Neighbouring voxel in the +z direction
    facetcountthisvoxel                    = facetcountthisvoxel+1;
    facetcountsofar                        = facetcountsofar+1;
    facetCOtemp(facetcountthisvoxel,1:3,1) = [ gridXupper(subX),gridYlower(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,2) = [ gridXupper(subX),gridYupper(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,3) = [ gridXlower(subX),gridYupper(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,4) = [ gridXlower(subX),gridYlower(subY),gridZupper(subZ) ];
    normalCOtemp(facetcountthisvoxel,1:3)  = [0,0,1];
  end
  if neighbourlist(loopP,4)==0   %Neighbouring voxel in the +y direction
    facetcountthisvoxel                    = facetcountthisvoxel+1;
    facetcountsofar                        = facetcountsofar+1;
    facetCOtemp(facetcountthisvoxel,1:3,1) = [ gridXupper(subX),gridYupper(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,2) = [ gridXlower(subX),gridYupper(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,3) = [ gridXlower(subX),gridYupper(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,4) = [ gridXupper(subX),gridYupper(subY),gridZupper(subZ) ];
    normalCOtemp(facetcountthisvoxel,1:3)  = [0,1,0];
  end
  if neighbourlist(loopP,5)==0   %Neighbouring voxel in the -z direction
    facetcountthisvoxel                    = facetcountthisvoxel+1;
    facetcountsofar                        = facetcountsofar+1;
    facetCOtemp(facetcountthisvoxel,1:3,1) = [ gridXlower(subX),gridYlower(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,2) = [ gridXlower(subX),gridYupper(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,3) = [ gridXupper(subX),gridYupper(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,4) = [ gridXupper(subX),gridYlower(subY),gridZlower(subZ) ];
    normalCOtemp(facetcountthisvoxel,1:3)  = [0,-1,0];
  end
  if neighbourlist(loopP,6)==0   %Neighbouring voxel in the +x direction
    facetcountthisvoxel                    = facetcountthisvoxel+1;
    facetcountsofar                        = facetcountsofar+1;
    facetCOtemp(facetcountthisvoxel,1:3,1) = [ gridXupper(subX),gridYupper(subY),gridZlower(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,2) = [ gridXupper(subX),gridYupper(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,3) = [ gridXupper(subX),gridYlower(subY),gridZupper(subZ) ];
    facetCOtemp(facetcountthisvoxel,1:3,4) = [ gridXupper(subX),gridYlower(subY),gridZlower(subZ) ];
    normalCOtemp(facetcountthisvoxel,1:3)  = [1,0,0];
  end
  
  facetsALL(facetcountsofar-facetcountthisvoxel+1:facetcountsofar,:,:) = facetCOtemp;
  normalsALL(facetcountsofar-facetcountthisvoxel+1:facetcountsofar,:)  = normalCOtemp;
  valuesALL(facetcountsofar-facetcountthisvoxel+1:facetcountsofar)     = gridINPUT(edgevoxelindices(loopP));

end

%======================================================
% PLOT THE QUADRANGULAR SURFACE MESH
%======================================================

if colYN==1
  
  %Define the colormap and resolution
  if exist('cmap','var')
    cres = size(cmap,1);
  else
    cres = 16;
    cmap = jet(cres);
  end

  if symYN==1
    %Ensure that the colorbar is symmetric about zero
    valueMINMAX = max(abs([min(valuesALL(:)),max(valuesALL(:))]));
    valuesTEMP  = valuesALL + valueMINMAX;
    valuesTEMP  = 1 + round( (cres-1)*valuesTEMP./(2*valueMINMAX) );    % The '1 + ...' ensures that the values are within [1:16] rather than [0:16]
    crange      = [-valueMINMAX,valueMINMAX];
  elseif symYN==0
    %The colorbar does not need to be symmetric about zero
    valuesTEMP  = valuesALL - min(valuesALL(:));
    valuesTEMP  = 1 + round( (cres-1)*valuesTEMP./max(valuesTEMP(:)) );    % The '1 + ...' ensures that the values are within [1:16] rather than [0:16]
    crange      = [min(valuesALL(:)),max(valuesALL(:))];
  end

  if crange(1) ~= crange(end)
    %If the facets have a range of colours
    hpat = cell(1,cres);
    for loopJ = 1:cres
      xco         = squeeze( facetsALL(valuesTEMP==loopJ,1,:) )';
      yco         = squeeze( facetsALL(valuesTEMP==loopJ,2,:) )';
      zco         = squeeze( facetsALL(valuesTEMP==loopJ,3,:) )';
      hpat{loopJ} = patch(xco,yco,zco,cmap(loopJ,:));
      set(hpat{loopJ},'EdgeColor',cmap(loopJ,:)/2);
    end
    set(gca, 'CLim', crange);
    if cbar==1
      hcbar = colorbar('NorthOutside','YLim',crange);
    elseif cbar==2
      hcbar = colorbar('EastOutside','YLim',crange);
    elseif cbar==3
      hcbar = colorbar('SouthOutside','YLim',crange);
    elseif cbar==4
      hcbar = colorbar('WestOutside','YLim',crange);
    end
    colormap(cmap);
  else
    %If the facets are all the same colour
    xco  = squeeze( facetsALL(:,1,:) )';
    yco  = squeeze( facetsALL(:,2,:) )';
    zco  = squeeze( facetsALL(:,3,:) )';
    hpat = patch(xco,yco,zco,cmap(round(end/2),:));
    set(hpat,'EdgeColor',cmap(round(end/2),:)/2);
    set(gca,'CLim',crange+[-1,1]);
    if cbar==1
      hcbar = colorbar('NorthOutside','YLim',crange+[-1,1]);
    elseif cbar==2
      hcbar = colorbar('EastOutside','YLim',crange+[-1,1]);
    elseif cbar==3
      hcbar = colorbar('SouthOutside','YLim',crange+[-1,1]);
    elseif cbar==4
      hcbar = colorbar('WestOutside','YLim',crange+[-1,1]);
    end
    colormap(cmap);
  end

elseif colYN==0
  
  %Define the colormap and resolution
  if exist('cmap','var')==1
    cmap = cmap(:)';  %Ensure colourmap is a 1x3 RGB triplet
  else
    cmap = [0,0,1];
  end

  %Plot a simple logical array where all facets are the same colour
  xco  = squeeze( facetsALL(:,1,:) )';
  yco  = squeeze( facetsALL(:,2,:) )';
  zco  = squeeze( facetsALL(:,3,:) )';
  hpat = patch(xco,yco,zco,cmap);
  set(hpat,'EdgeColor',cmap/2);

end
axis equal tight;


%======================================================
% SET THE OUTPUT PARAMETERS
%======================================================

if nargout>=1
  varargout(1) = {hpat};
end
if nargout==2
  if exist('hcbar','var')==1
    varargout(2) = {hcbar};
  else
    varargout(2) = {[]};
  end
end
