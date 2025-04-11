function [ind_transform, indgridlims] = transformToWindows(tdata, tgrid, windowlen)
%% [ind_transform, indgridlims] = TRANSFORMTOWINDOWS(tdata, tgrid, windowlen)
%
%   inputs
%       - tdata: vector with (gridded) independent variable (e.g. time) of data.
%       - tgrid: euqally spaced (e.g. time) grid for windowed calculations.
%       - windowlen: length of window.
%
%   outputs
%       - ind_transform: a 2xN array with first and last indices of each of
%                        the N windows of tdata that are within tgrid.
%       - indgridlims: 1x2 double first and last indices of tgrid that
%                      fully contain tdata in their windows.
%
%
% TRANSFORMTOWINDOWS.m looks at two grids (e.g. time vectors) and returns
% the information to reshape a vector into a matrix. Then, column-wise
% operations can be calculated efficiently on the data that is reshaped
% into a matrix.
%
% TRANSFORMTOWINDOWS.m is useful for timeseries analysis, such as computing
% a moving average, computing multiple spectra, etc. TRANSFORMTOWINDOWS.m
% is similar to reshapeintoWindows.m, but here windowlen is independent of
% the spacing in tgrid, and can include overlap independent 
%
% TRANSFORMTOWINDOWS.m requires that (a) windowlen is a multiple of data
% sampling rate and (b) there must data points that match with grid points.
%
% See also:
%   RESHAPEINTOWINDOWS.m
%
% Olavo Badaro Marques.


%% Check tdata and tgrid are both equally spaced grids

%
dt_data = diff(tdata(1:2));
% % %
% % if any(diff(tdata) ~= dt_data)
% %     error('Data grid is not an equally spaced grid.')
% % end

% % %
% % windowlen = diff(tgrid(1:2));
% % %
% % if any(diff(tgrid) ~= windowlen)
% %     error('Data grid is not an equally spaced grid.')
% % end


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


%%

%
Nfullsegments = indgridlims(2) - indgridlims(1) + 1;

%
if mod(npts_window, 2)==0
    %
    vec_inds = (-npts_window/2) : 1 : (npts_window/2 - 1);
else
    %
    vec_inds = ((-npts_window/2) + 0.5) : 1 : ((npts_window/2) - 0.5);
end

%
vec_inds_bounds = vec_inds([1, end]);

%%

%
ind_tdata_ongrid_first = find(tdata == tgrid(indgridlims(1)));
ind_tdata_ongrid_last = find(tdata == tgrid(indgridlims(2)));

%
if indgridlims(1) == indgridlims(2)
    
    %
    ind_dataongrid = ind_tdata_ongrid_first;
else
    
    %
    ind_tdata_ongrid_second = find(tdata == tgrid(indgridlims(1)+1));

    %
    ind_diff = ind_tdata_ongrid_second - ind_tdata_ongrid_first;

    %
    ind_dataongrid = ind_tdata_ongrid_first : ind_diff : ind_tdata_ongrid_last;
end




%%

%
if ind_dataongrid(end)~=ind_tdata_ongrid_last
    warning('error!!!')
    keyboard
end

%
ind_transform = repmat(ind_dataongrid(:).', 2, 1) + ...
                repmat(vec_inds_bounds(:), 1, Nfullsegments);

