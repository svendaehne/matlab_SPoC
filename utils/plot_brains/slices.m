function[x] = slices(m,d)
%SLICES  return a cell array of slices of a tensor along the specied dimension
%
% This function is useful for slicing high-dimensional images (e.g. fMRI
% images) into lower-dimensional pieces.  The output can then be either
% processed directly, or passed into cellfun or arrayfun for further
% processing.
%
% USAGE:
%   x = slices(m,[d])
%
% INPUTS:
%
%   m: a tensor (i.e. multi-dimensional matrix) that you wish to slice
%
%   d: the dimnsion along which you wish to make the slices.  default: 
%      d = ndims(m).
%
% OUTPUTS:
%
%   x: a 1 by size(m,d) cell array of slices of m.  each element of the
%      cell array is an  (ndims(m) - 1)-dimensional tensor.
%
% EXAMPLE:
%
%   %create a tensor to slice
%   m = reshape(1:100,10,10);
%
%   %slice along first dimension: returns rows of m
%   s1 = slices(m,1);
%   
%   %slice along second dimension: returns columns of m
%   s2 = slices(m,2);
%
%   %slice along third (or greater) dimension: returns m as a cell array
%   s3 = slices(m,3);   
%
% SEE ALSO: REPMAT, CELLFUN, ARRAYFUN, JOIN
%
%   AUTHOR: JEREMY R. MANNING
%  CONTACT: manning3@princeton.edu


%CHANGELOG
%3-3-12    JRM      Wrote it.

%returns slices along the d^th (or last) dimension of m
if ~exist('d','var'), d = ndims(m); end
x = arrayfun(@slice_matrix,repmat({m},1,size(m,d)),d*ones(1,size(m,d)),1:size(m,d),'UniformOutput',false);



%return the ith slice along the d^th dimension of m
function[x] = slice_matrix(m,d,i)
if iscell(m) && (length(m) == 1)
    m = m{1};
end
x = eval(sprintf(['m(',indexer(d,ndims(m)),')'],i));

%i: index to change; n: number of dimensions in m
function[s] = indexer(i,n) 
n = max([i n]);
if max([n i]) < 1
    s = '';
    return;
end
%make a string of n ":"'s, separated by commas
s = join(',',repmat({':'},1,n));

%replace the ith ':' with '%d'
inds = strfind(s,':');
sep = inds(i);
s = [s(1:sep-1),'%d',s(sep+1:end)];


%COPIED FROM JOIN BY GERALD DALLEY: http://www.mathworks.com/matlabcentral/fileexchange/4872-join
function s = join(d,varargin)
%S=JOIN(D,L) joins a cell array of strings L by inserting string D in
%            between each element of L.  Meant to work roughly like the
%            PERL join function (but without any fancy regular expression
%            support).  L may be any recursive combination of a list 
%            of strings and a cell array of lists.
%
%For any of the following examples,
%    >> join('_', {'this', 'is', 'a', 'string'} )
%    >> join('_', 'this', 'is', 'a', 'string' )
%    >> join('_', {'this', 'is'}, 'a', 'string' )
%    >> join('_', {{'this', 'is'}, 'a'}, 'string' )
%    >> join('_', 'this', {'is', 'a', 'string'} )
%the result is:
%    ans = 
%        'this_is_a_string'
%
%Written by Gerald Dalley (dalleyg@mit.edu), 2004
if isempty(varargin)
    s = '';
else
    if (iscell(varargin{1}))
        s = join(d, varargin{1}{:});
    else
        s = varargin{1};
    end
    
    for ss = 2:length(varargin)
        s = [s d join(d, varargin{ss})]; %#ok<AGROW>
    end
end
