function [Sxy, frequency, indspeclims, dof] = compute_cpsd(x, y, dtime_data, dtime_spec, dt_window, dt_fft, fracoverlap, taperwindow)
%% [Sxy, frequency, indspeclims, dof] = COMPUTE_CPSD(x, y, dtime_data, dtime_spec, dt_window, dt_fft, fracoverlap, taperwindow)
%
%   inputs
%       - x: one timeseries (with evenly spaced data).
%       - y: a second timeseries.
%       - dtime_data: datetime of the data x and y.
%       - dtime_spec: datetime for computing the cross-spectra.
%       - dt_window: duration variable for ensemble-averaging spectra
%                    (i.e., the window around dtime_spec grid points).
%       - dt_fft: duration variable with the time length for computing
%                 fft to data segments.
%       - fracoverlap: a scalar >=0 and <1, that gives the fraction of
%                      overlap between data segments that are fft'ed.
%       - taperwindow:
%
%   outputs
%       - Sxy: cross-spectra between x and y.
%       - frequency: frequency vector.
%       - indspeclims: indices that give the first and last indices of
%                      dtime_spec, where Sxy could be computed.
%       - dof: vector with degrees of freedom for each spectrum in sxy.
%
%
% COMPUTE_CPSD.m is a high-level function for computing cross-spectra Sxy
% between x and y. The cross-spectra are coputed by Matlab's cpsd.m. Most
% of COMPUTE_CPSD.m is dedicated to subsetting the data (x and y) for
% computing Sxy at the corresponding dtime_spec, given the other parameters
% (e.g., dt_window, dt_fft, etc).
%
% To subset the data, either the function reshapeintoWindows.m or
% transformToWindows.m is used.  
%
%
% See also:
%   cpsd.m
%   window.m
%   reshapeintoWindows.m
%   transformToWindows.m


%% Check if there are any data points within the limits of dtime_spec.
% 

%
if length(dtime_spec)>1
    if ~any((dtime_data >= dtime_spec(1)) & (dtime_data <= dtime_spec(end)))
        %
        warning(['Error: no data available for computing ' ...
                 'spectra at datetimes dtime_spec'])
        %
        Sxy = NaN;
        frequency = NaN;
        indspeclims = NaN;
        %
        return
    end
end


%% Look at dtime_spec and dt_window to see how
% data vectors will be changed to matrices

%
if length(dtime_spec) > 1
    %
    if all(diff(dtime_spec)==dt_window)
        lreshape = true;
    else
        lreshape = false;
    end

%
elseif isscalar(dtime_spec)
    %
    lreshape = false;
else
    %
    error('unexpected dtime_spec')
    
end


%% Find indices to rearrange data so that data
% for each dtime_spec(i) is in a column vector

%
if lreshape
    %
    [indsub, reshapeNdims, indgridlims] = reshapeintoWindows(dtime_data, dtime_spec);
    
else
    %
    if length(dtime_spec) > 1
        % 
        [ind_transform, indgridlims] = transformToWindows(dtime_data, dtime_spec, dt_window);

	%
    else   

        %
        ind_1 = find(dtime_data >= (dtime_spec - (dt_window/2)), 1, 'first');
        ind_2 = find(dtime_data  < (dtime_spec + (dt_window/2)), 1, 'last');
        %
        ind_transform = [ind_1; ind_2];
        indgridlims = 1;
    end
    
end

%
indspeclims = indgridlims;


%% Shift data vectors x and y to matrices

% If possible, only reshape
if lreshape
    %
    x_array = reshape(x(indsub), reshapeNdims);
    y_array = reshape(y(indsub), reshapeNdims);

    %
    N_ptstime_sub = reshapeNdims(2);
    
% Or transform vectors if it needs to be more complicated
else
    %
    N_cols = size(ind_transform, 2);
    N_rows = ind_transform(2, 1) - ind_transform(1, 1) + 1;
    %
    x_array = NaN(N_rows, N_cols);
    y_array = x_array;
    %
    for i = 1:N_cols
        
        %
        inds_get = ind_transform(1, i) : 1 : ind_transform(2, i);
        %
        x_array(:, i) = x(inds_get);
        y_array(:, i) = y(inds_get);
    end
    
    %
    N_ptstime_sub = N_cols;
end


%% Get number of points per FFT segment

%
dt_data = seconds(diff(dtime_data(1:2)));

%
if ~isduration(dt_fft)
    error('wrong input format!')
end

%
nfft = seconds(dt_fft) * (1/dt_data);

%
freqNyquist = (1/dt_data)/2;
df = 1/(nfft*dt_data);
%
Nfreq = length(0 : df : freqNyquist);


%% Loop over columns of x_array (and y_array), reshape so that each
% column is a segment that will be fft'ed, then detrend each segment and
% apply the fft
 
%
nstep_chunk = floor(nfft*(1-fracoverlap));
%
ind_first_chunk = 1:nstep_chunk:size(x_array, 1);
ind_last_chunk = ind_first_chunk + nfft - 1;
%
linsubsegment = (ind_last_chunk <= size(x_array, 1));
%
ind_first_chunk = ind_first_chunk(linsubsegment);
ind_last_chunk = ind_last_chunk(linsubsegment);
%
Nchunks = length(ind_last_chunk);

%
dof_max = 2*Nchunks;
%
dof_TH = floor(0.8*dof_max);

%
ind_get_chunks = cell(1, Nchunks);
%
for i = 1:Nchunks
    ind_get_chunks{i} = ind_first_chunk(i) : 1 : ind_last_chunk(i);
end

%
Sxy_new = NaN(N_ptstime_sub, Nfreq);
dof = zeros(N_ptstime_sub, 1);

%
for i1 = 1:N_ptstime_sub
    
    %%

    %
    x_chunks = NaN(nfft, Nchunks);
    y_chunks = x_chunks;
    
    %
    for i2 = 1:Nchunks
        %
        x_chunks(:, i2) = x_array(ind_get_chunks{i2}, i1);
        y_chunks(:, i2) = y_array(ind_get_chunks{i2}, i1);
    end
    
    %
    lok_chunk = all(~isnan(x_chunks), 1);
    %
    dof_aux = 2*length(lok_chunk(lok_chunk));
    
    %
    if dof_aux < dof_TH
        continue
    end
    
    
    %%
    
    % Detrend
    x_chunks = detrend(x_chunks, 1, 'omitnan');
    y_chunks = detrend(y_chunks, 1, 'omitnan');
    
    
    % Now compute cpsd
    [Sxy_chunks, frequency] = cpsd(x_chunks(:, lok_chunk), y_chunks(:, lok_chunk), ...
                                   window(taperwindow, nfft), ...
                                   0, ...
                                   nfft, 1./dt_data);
    
	%
    Sxy_new(i1, :) = mean(Sxy_chunks, 2);

    %
	dof(i1) = dof_aux;


    
end

%
Sxy = Sxy_new;

