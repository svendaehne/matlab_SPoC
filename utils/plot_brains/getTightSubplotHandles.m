function[hs] = getTightSubplotHandles(FIG_MARGIN,ROW_MARGIN,COLUMN_MARGIN,N_ROWS,N_COLUMNS)
%GETTIGHTSUBPLOTHANDLES  return handles to subplots of with the specified spacing
%
% syntax: hs = getTightSubplotHandles(outer_margin,row_margin,...
%                                     column_margin,n_rows,n_columns);
%
% INPUTS:
%
%   outer_margin: specifies proportion (between 0 and 1) of the figure to
%                 use as a border.  outer_margin = 0 uses the full figure
%                 area (no border).
%
%     row_margin: specifies proportion (between 0 and 1) of the figure to
%                 use as a border between subplot rows.  row_margin = 0
%                 sets rows to be directly adjacent (no border).
%
%  column_margin: specifies proportion (between 0 and 1) of the figure to
%                 use as a border between subplot columns.  column_margin =
%                 0 sets columns to be directly adjacent (no border).
%
%         n_rows: number of subplot rows
%
%      n_columns: number of subplot columns
%
% OUTPUTS:
%
%             hs: handles to the subplots (a 1 by n_rows*n_columns vector)
%
% EXAMPLE:
%
%   %create subplots with 1% margins between rows, 2% margins between
%   %columns, and a 5% border. create 15 subplots, arranged in 5 rows and 3
%   %columns.
%   hs = getTightSubplotHandles(0.05,0.01,0.02,5,3);
%
%   %plot something in each subplot
%   for i = 1:length(hs)
%       axes(hs(i));
%       plot(rand(1,10),'k','LineWidth',2);
%       axis off;
%   end
%
% SEE ALSO: PLOT, SUBPLOT, AXES, AXIS
%
%   AUTHOR: JEREMY R. MANNING
%  CONTACT: manning3@princeton.edu

%CHANGELOG
% 3-3-12    JRM    Wrote it.

assert((FIG_MARGIN >= 0) & (FIG_MARGIN < 1),'specify outer margin between 0 and 1');
assert((ROW_MARGIN >= 0) & (ROW_MARGIN < 1),'specify row margin between 0 and 1');
assert((COLUMN_MARGIN >= 0) & (COLUMN_MARGIN < 1),'specify row margin between 0 and 1');
assert(N_ROWS > 0,'number of rows must be positive');
assert(N_COLUMNS > 0,'number of columns must be positive');

w = 1; %full figure width
h = 1; %full figure height

w = w-(2*FIG_MARGIN)-((N_COLUMNS-1)*COLUMN_MARGIN); %workable width
h = h-(2*FIG_MARGIN)-((N_ROWS-1)*ROW_MARGIN); %workable height

sw = w/N_COLUMNS; %subplot width
sh = h/N_ROWS; %subplot height

hs = zeros(1,N_COLUMNS*N_ROWS);
clf;
count = 1;
startH = 1 - FIG_MARGIN - sh;
for i = 1:N_ROWS    
    startW = FIG_MARGIN;
    for j = 1:N_COLUMNS        
        h = subplot('position',[startW startH sw sh]);
        hs(count) = h;
        startW = startW + sw + COLUMN_MARGIN;        
        count = count+1;
    end
    startH = startH - sh - ROW_MARGIN;
end
