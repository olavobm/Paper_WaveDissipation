function [a1, b1, a2, b2, frequency, Szz, indgridlims] = Spotter_proc_dirmoments(spotterL1, dtime_grid, fftparams)
%% [a1, b1, a2, b2, frequency, Szz, indgridlims] = SPOTTER_PROC_DIRMOMENTS(spotterL1, dtime_grid, fftparams)
%
%   inputs
%       - spotterL1: L1 spotter data structure.
%       - dtime_grid: time grid in datetime.
%       - fftparams: structure variable where fields have the information
%                    to compute spectra.
%
%   outputs
%       - a1, b1: first directional moments.
%       - a2, b2: second directional moments.
%       - frequency: frequency vector (in Hz).
%       - Szz: frequency spectra of the vertical sea-surface elevation.
%       - indgridlims: 1x2 array with first and last indices of dtime_grid
%                      that fully contain data for spectra.
%
%
% SPOTTER_PROC_DIRMOMENTS.m is a high-level function that computes
% directional wave moments from level 1 Spotter data in the variable
% structure spotterL1.
%
%
% See also:
%   wave_dirmoments_xyz.m


%% Pass parameters to shorter variables

%
dt_window = fftparams.dt_avgwindow;
dt_fft = fftparams.dt_fft;
fftoverlap = fftparams.fftoverlap;
ffttaper = fftparams.ffttaper;


%% Compute directional moments

%
[frequency, ...
 a1, a2, b1, b2, ...
 Szz, indgridlims] = wave_dirmoments_xyz(spotterL1.displacement.x, ...
                                         spotterL1.displacement.y, ...
                                         spotterL1.displacement.z, ...
                                         spotterL1.displacement.dtime, ...
                                         dtime_grid, ...
                                         dt_window, dt_fft, ...
                                         fftoverlap, ffttaper);

                                     
%% NaN-out results below low frequency cut-off

%
lbelowcutoff = (frequency <= fftparams.lowfreq_cutoff);
%
a1(:, lbelowcutoff) = NaN;
a2(:, lbelowcutoff) = NaN;
b1(:, lbelowcutoff) = NaN;
b2(:, lbelowcutoff) = NaN;
Szz(:, lbelowcutoff) = NaN;

