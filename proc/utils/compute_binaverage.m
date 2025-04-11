function [xmean, indgridlims] = compute_binaverage(x, dtime_data, dtime_avg, dt_window)
%% [xmean, indgridlims] = COMPUTE_BINAVERAGE(x, dtime_data, dtime_avg, dt_window)
%
%   inputs
%       - x: data vector.
%       - dtime_data: datetime vector corresponding to x.
%       - dtime_avg: datetime vector for computing means.
%       - dt_window: duration variable for window length to compute the mean.
%
%   outputs
%       - xmean: bin-averaged x.
%       - indgridlims: indices that relate xmean to the corresponding grid
%                      points of dtime_avg.
%
%
% COMPUTE_BINAVERAGE.m computes bin averages of x. Most of what this
% function does is to get the appropriate data for computing means
% correspondent to times dtime_avg. This is done by function
% reshapeintoWindows.m or transformToWindows.m.
%
%
% See also:
%   reshapeintoWindows.m
%   transformToWindows.m


%% Look at dtime_avg and dt_window to see how
% data vectors will be changed to matrices

%
if length(dtime_avg) > 1
    %
    if all(diff(dtime_avg)==dt_window)
        lreshape = true;
    else
        lreshape = false;
    end

%
elseif isscalar(dtime_avg)
    %
    lreshape = false;
else
    %
    error('unexpected dtime_avg')
    
end


%% Get indices for turning vector into a matrix

%
if lreshape
    %
    [indsub, reshapeNdims, indgridlims] = reshapeintoWindows(dtime_data, dtime_avg);
    
else
    %
    [ind_transform, indgridlims] = transformToWindows(dtime_data, dtime_avg, dt_window);
    
end


%% Shift data vector x to a mtrix

% If possible, only reshape
if lreshape
    %
    x_array = reshape(x(indsub), reshapeNdims);
    
% Or transform vectors if it needs to be more complicated
else
    %
    N_cols = size(ind_transform, 2);
    N_rows = ind_transform(2, 1) - ind_transform(1, 1) + 1;
    %
    x_array = NaN(N_rows, N_cols);
    %
    for i = 1:N_cols
        %
        inds_get = ind_transform(1, i) : 1 : ind_transform(2, i);
        %
        x_array(:, i) = x(inds_get);
    end
    
end


%% Compute (column-wise) mean

%
xmean = mean(x_array, 1);


