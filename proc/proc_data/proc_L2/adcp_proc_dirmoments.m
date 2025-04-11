function [a1, b1, a2, b2, frequency, Szz, indgridlims, dof_p, dof_u] = adcp_proc_dirmoments(dataL1, dtime_grid, fftparams)
%% [a1, b1, a2, b2, frequency, Szz, indgridlims, dof_p, dof_u] = ADCP_PROC_DIRMOMENTS(dataL1, dtime_grid, fftparams)
%
%   inputs
%       - dataL1: structure variable with ADCP data.
%       - dtime_grid: vector datetime grid for computing spectra.
%       - fftparams: structure variable with parameters for computing
%                    spectra.
%
%   outputs
%       - a1, b1: first directional moments.
%       - a2, b2: second directional moments.
%       - frequency: wave frequency in Hz.
%       - Szz: frequency spectra of pressure.
%       - indgridlims: 1x2 array with first and last indices of dtime_grid
%                      that fully contain data for spectra.
%       - dof_p: vector with degrees of freedom of the pressure spectra.
%       - dof_u: vector with degrees of freedom of the velocity spectra.
%
%
% ADCP_PROC_DIRMOMENTS.m is a high-level function that computes directional
% wave moments from level 1 ADCP data in the variable structure dataL1.
%
%
% See also:
%   wave_dirmoments_uvp.m


%% Pass parameters to shorter variables

dt_window = fftparams.dt_avgwindow;
dt_fft = fftparams.dt_fft;
fftoverlap = fftparams.fftoverlap;
ffttaper = fftparams.ffttaper;


%% Compute directional moments

%
[frequency, ...
 a1, a2, b1, b2, ...
 Szz, indgridlims, ...
 dof_p, dof_u] = wave_dirmoments_uvp(dataL1.u(:), ...
                                     dataL1.v(:), ...
                                     dataL1.pressure, ...
                                     dataL1.dtime, ...
                                     dtime_grid, ...
                                     dt_window, dt_fft, fftoverlap, ffttaper);

                                     
% definitions are such that atan2(b1, a1) is direction where waves
% propagate to in terms of the counterclockwise angle from the u>0
% direction (i.e. the standard trigonometric convention).

                                     
%% NaN-out low frequencies that are unreliable because of cpsd

%
lbelowcutoff = (frequency <= fftparams.lowfreq_cutoff);
%
a1(:, lbelowcutoff) = NaN;
a2(:, lbelowcutoff) = NaN;
b1(:, lbelowcutoff) = NaN;
b2(:, lbelowcutoff) = NaN;
Szz(:, lbelowcutoff) = NaN;


%% Trim out high-frequencies

%
lbelow_highcutoff = (frequency <= fftparams.highfreq_cutoff);

%
frequency = frequency(lbelow_highcutoff);
a1 = a1(:, lbelow_highcutoff);
b1 = b1(:, lbelow_highcutoff);
a2 = a2(:, lbelow_highcutoff);
b2 = b2(:, lbelow_highcutoff);
Szz = Szz(:, lbelow_highcutoff);




