function [frequency, a1, a2, b1, b2, Szz, indgridlims] = wave_dirmoments_xyz(x, y, z, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper)
%% [frequency, a1, a2, b1, b2, Szz, indgridlims] = WAVE_DIRMOMENTS_XYZ(x, y, z, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper)
%
%   inputs
%       - x, y, z: horizontal (x, y) and vertical (z) displamecenets.
%       - dtime_data: gridded datetime of the displacement data.
%       - dtime_grid: time grid (in datetime) to calculate dir. moments.
%       - dt_window: time window for each spectral estimate.
%       - dt_fft: time window for the individual fft segments.
%       - fftoverlap: overlap fraction of the fft segments.
%       - ffttaper: Matlab handle of the window to taper fft segments.
%
%   outputs
%       - frequency: frequency vector (in Hz).
%       - a1, a2, b1, b2: directional moments.
%       - Szz: sea-surface elevation spectrum (a).
%       - indgridlims: 1x2 array with first and last indices of dtime_grid
%                      that that fully contain dtime_data in their windows.
%
%
% WAVE_DIRMOMENTS_XYZ.m takes displacement data (e.g. from Spotter
% wave buoy) and calculates directional moments a1, b1, a2, and b2.
% Spectra (and cross-spectra) are calculated through Matlab's cpsd.m
% function, but through the high-level function compute_cpsd.m
%
% See also:
%   COMPUTE_CPSD.m


%% Compute spectra and cross-spectra

%
[Szz, ...
 frequency, ...
 indgridlims] = compute_cpsd(z, z, dtime_data, dtime_grid, ...
                                   dt_window, dt_fft, fftoverlap, ffttaper);

%
Sxx = compute_cpsd(x, x, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
%
Syy = compute_cpsd(y, y, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
%
Sxz = compute_cpsd(x, z, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
%
Syz = compute_cpsd(y, z, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
%
Sxy = compute_cpsd(x, y, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);


%% Compute direction moments from wave buoy displacements

% This is such that atan2(b1, a1) gives the direction
% of where waves are propagating to
a1 = -imag(Sxz) ./ (sqrt(Szz.*(Sxx + Syy)));
b1 = -imag(Syz) ./ (sqrt(Szz.*(Sxx + Syy)));
%
a2 = (Sxx - Syy)./(Sxx + Syy);
b2 = (2.*real(Sxy)) ./ (Sxx + Syy);


