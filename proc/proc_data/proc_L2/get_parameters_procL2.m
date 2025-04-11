function [dtime_grid, fftparams] = get_parameters_procL2()
%% [dtime_grid, fftparams] = GET_PARAMETERS_PROCL2()
%
%   outputs
%       - dtime_grid: datetime vector.
%       - fftparams: structure variable with parameters
%                    for computing spectra.
%
%
% GET_PARAMETERS_PROCL2.m simply returns the dtime grid vector and the
% parameters used to compute spectra in the same manner for all instruments
% in the L2 data processing.


%% Define time-grid for computing spectra

%
dtime_lims = [datetime(2022, 06, 15, 12, 00, 00), ...
              datetime(2022, 07, 21, 05, 00, 00)];
%
dtime_lims.TimeZone = 'America/Los_Angeles';

%
dt_grid = hours(1);
dtime_grid = dtime_lims(1) : dt_grid : dtime_lims(2);
%
dtime_grid = dtime_grid(:);


%% Define parameters for spectra

% FFT window (in seconds)
dt_fft = 120;

% Averaging window (in seconds)
dt_avgwindow = seconds(dt_grid);

%
fft_overlap = 0.5;

%
fft_taperwindow = @hann;


%% Define a frequency cutoff to NaN-out spectral estimates
% lower than a certain (low) frequency

%
lowfreq_cutoff = 0.04;


%% Define a high frequency cutoff to remove spectral estimates
% at very high frequencies because they may not be reliable
% (e.g. large instrument noise).
%
% Note that these high frequencies will be trimmed out of the
% spectral arrays, instead of just NaNed-out.

%
highfreq_cutoff = 1.25;    % 1.25 Hz is the Spotter Nyquist frequency


%% Put the parameters above in a structure

%
fftparams.dt_fft = seconds(dt_fft);    % in duration (seconds) format
fftparams.fftoverlap = fft_overlap;
fftparams.ffttaper = fft_taperwindow;
fftparams.dt_avgwindow = seconds(dt_avgwindow);    % in duration (seconds) format
%
fftparams.lowfreq_cutoff = lowfreq_cutoff;
fftparams.highfreq_cutoff = highfreq_cutoff;

