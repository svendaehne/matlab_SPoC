function p = sanePColor(varargin)
%SANEPCOLOR  simple wrapper for pcolor
%
% Unlike the built-in pcolor command, this function does not cut off the
% last row and column of the input matrix.  In this way, sanePColor is
% intended to be as easy to use as imagesc, but allows the user to specify
% the x and y coordinates of each cell if desired.  This function is also
% useful as an alternative means of generating images to print to PDF that
% are compatible with OS X's "Preview" PDF viewer (imagesc images appear
% "blurred" when printing to a PDF as a vector graphic and viewed using
% Preview).
%
% Usage: p = sanePColor([x,y],z,[logx],[logy]);
%
%INPUTS:
%
%    x: an array of x values.  can also specify a min and max x value.
%       these values correspond to columns of z. [IF THIS ARGUMENT IS USED,
%       MUST ALSO SPECIFY Y VALUES.]
% 
%    y: an array of y values.  can also specify a min and max y value.
%       these values correspond to rows of z.  [IF THIS ARGUMENT IS USED,
%       MUST ALSO SPECIFY X VALUES.]
% 
%    z: a 2d matrix of values.  this matrix determines the color at each
%       point.
% 
% logx: if this optional argument is set to true, the x-axis will plotted
%       in log scale (similar to semilogx).
%
% logy: if this optional argument is set to true, the y-axis will plotted
%       in log scale (similar to semilogy).
%
%OUTPUTS:
%
%    p: a handle to the resulting pcolor image.
%
% EXAMPLE:
%
%   m = membrane;
%   p = sanePColor(m);
%
% SEE ALSO: PCOLOR, IMAGE, IMAGESC, SEMILOGX, SEMILOGY, LOGLOG, PADARRAY
%
%   AUTHOR: JEREMY R. MANNING
%  CONTACT: manning3@princeton.edu


%CHANGELOG
%3-16-10    JRM      Wrote it.
%3-12-12    JRM      Support a more diverse range of input configurations.

%parse arguments
if length(varargin) == 1 %just a z
    z = varargin{1};
    x = [1 size(z,2)];
    y = [1 size(z,1)];
    [logx,logy] = deal(false);    
elseif (length(varargin) >= 4) %x, y, z, and possibly logx and/or logy
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    logx = varargin{4};
    if length(varargin) >= 5, logy = varargin{5}; else logy = false; end
elseif length(varargin) == 2 %z and logx
    z = varargin{1};
    logx = varargin{2};
    logy = false;
    x = [1 size(z,2)];
    y = [1 size(z,1)];
else %length(varargin) == 3
    if isempty(varargin)
        fprintf('\nUsage: p = sanePColor([x,y],z,[logx],[logy]);\n');
        fprintf('Type ''help %s'' for more info.\n\n',mfilename);
        p = [];
        return;
    end
    %posibility 1: x, y, z
    if length(varargin{1}) > 1 && length(varargin{2}) > 1
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        [logx,logy] = deal(false);
    %posibility 2: z, logx, and logy
    else
        z = varargin{1};
        x = [1 size(z,2)];
        y = [1 size(z,1)];
        logx = varargin{2};
        logy = varargin{3};
    end
end

z = padarray(z,[1 1],'replicate','post');
if logx
    newx = logspace(log10(min(x)),log10(max(x)),size(z,2));
else
    newx = linspace(min(x),max(x),size(z,2));
end

if logy
    newy = logspace(log10(min(y)),log10(max(y)),size(z,1));
else
    newy = linspace(min(y),max(y),size(z,1));
end

p = pcolor(newx,newy,z);
set(gca,'XTickMode','Manual','YTickMode','Manual');
if exist('logx','var') && logx
    set(gca,'XScale','log');
    set(gca,'XTick',logspace(log10(min(x)),log10(max(x)),4));
else
    set(gca,'XTick',linspace(min(x),max(x),4));
end
if exist('logy','var') && logy
    set(gca,'YScale','log');
    set(gca,'YTick',logspace(log10(min(y)),log10(max(y)),8));
else
    set(gca,'YTick',linspace(min(y),max(y),8));
end

shading flat;