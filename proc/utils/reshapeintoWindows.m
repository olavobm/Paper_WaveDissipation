function [indsub, reshapeNdims, indgridlims] = reshapeintoWindows(tdata, tgrid)
%% [indsub, reshapeNdims, indgridlims] = RESHAPEINTOWINDOWS(tdata, tgrid)
%
%   inputs
%       - tdata: vector with (gridded) independent variable (e.g. time) of data.
%       - tgrid: equally spaced (e.g. time) grid for windowed calculations
%                (where the window length is diff(tgrid)).
%
%   outputs
%       - indsub: indices of tdata that are between the windows of the
%                 first and last grid points tgrid([1, end]).
%       - reshapeNdims: 1x2 double with rowsxcolumns to reshape
%                       vector data into a a matrix.
%       - indgridlims: 1x2 double first and last indices of tgrid that
%                      fully contain tdata in their windows.
%
%
% RESHAPEINTOWINDOWS.m looks at two grids (e.g. time vectors) and returns
% the information to reshape a vector into a matrix, where each column of
% the latter corresponds to a grid point of tgrid. Then column-wise
% operations can be calculated efficiently on the data that is reshaped
% into a matrix.
%
% RESHAPEINTOWINDOWS.m is useful for timeseries analysis, such as computing
% a moving average, computing multiple spectra, etc.
%
% RESHAPEINTOWINDOWS.m is only defined for a grid tgrid where the windows
% of the operation are equal to the grid spacing (i.e. no overlap of the
% windows is allowed). Because of this, tgrid must have at least length==2.
%
%
% Olavo Badaro Marques.


%% Check tdata and tgrid are both equally spaced grids

%
dt_data = diff(tdata(1:2));
%
if any(diff(tdata) ~= dt_data)
    error('Data grid is not an equally spaced grid.')
end

%
windowlen = diff(tgrid(1:2));
%
if any(diff(tgrid) ~= windowlen)
    error('Data grid is not an equally spaced grid.')
end


%% Check tdata and tgrid have the same class

if ~isequal(class(tdata), class(tgrid))
    error('!!!!')
end


%% Get the number of data points per window of tgrid

%
if isdatetime(tdata)
    %
    npts_window = (1/seconds(dt_data)) * seconds(windowlen);
%
else
    %
    npts_window = (1/dt_data) * (windowlen);
end


%% Get the first and last grid points that
% have windows entirely covered by tdata

%
tgrid_bounds = [(tgrid(:) - (windowlen/2)).'; ...
                (tgrid(:) + (windowlen/2)).'];

%
ind_first_wholeinterval = find(tgrid_bounds(1, :) >= tdata(1), 1, 'first');
ind_last_wholeinterval = find(tgrid_bounds(2, :) < tdata(end), 1, 'last');

%
indgridlims = [ind_first_wholeinterval, ind_last_wholeinterval];


%% Get indices of tdata that correspond to tgrid(indgridlims)

%
ind_data_start_firstwindow = find(tdata == tgrid_bounds(1, ind_first_wholeinterval));

%
ind_data_end_lastwindow = find(tdata == tgrid_bounds(2, ind_last_wholeinterval));
ind_data_end_lastwindow = ind_data_end_lastwindow - 1;

%
indsub = ind_data_start_firstwindow : ind_data_end_lastwindow;


%% Get the dimensions to reshape the vector tdata to a matrix where
% each column corresponds to a grid point of tgrid

%
Nindsub = length(indsub);

% This should not happen. Still it's checking for an error
if mod(Nindsub, npts_window)~=0
    error('Integer not obtained. Reshape failed.')
end

%
npts_intervals = Nindsub/npts_window;

%
reshapeNdims = [npts_window, npts_intervals];


