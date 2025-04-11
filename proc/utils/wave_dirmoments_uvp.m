function [frequency, a1, a2, b1, b2, Spp, indgridlims, dof_p, dof_u] = wave_dirmoments_uvp(u, v, p, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper)
%% [frequency, a1, a2, b1, b2, Spp, indgridlims, dof_p, dof_u] = WAVE_DIRMOMENTS_UVP(u, v, p, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper)
%
%   inputs
%       - u: horizontal velocity component in the x direction.
%       - v: horizontal velocity component in the y direction.
%       - p: pressure.
%       - dtime_data: datetime of the data (u, v, p).
%       - dtime_grid: datetime grid of where to compute a1, b1, a2, and b2.
%       - dt_window: window length for averaging statistics. 
%       - dt_fft: window length for computing FFT.
%       - fftoverlap: scalar (< 1) with the overlap between FFT windows.
%       - ffttaper: window handle for tapering data in FFT windows.
%
%   outputs
%       - frequency: frequency in Hz.
%       - a1, a2, b1, b2: directional moments.
%       - Spp: pressure spectra.
%       - indgridlims: 1x2 array with first and last indices of dtime_grid
%                      that fully contain dtime_data in their windows.
%       - dof_p: degrees of freedom in pressure spectra.
%       - dof_u: degrees of freedom in velocity spectra.
%
%
% WAVE_DIRMOMENTS_UVP.m takes ADCP data and calculates first and second
% directional moments (a1, b1, a2, and b2).
%
% Cross-spectra is calculated by Matlab's cpsd.m, but through my high-level
% function compute_cpsd.m
%
% See also:
%   COMPUTE_CPSD.m
%
% Written by Olavo Marques Oct/2022.


%% Compute cross-spectra

%
[Spp, ...
 frequency, ...
 indgridlims, ...
 dof_p] = compute_cpsd(p, p, dtime_data, dtime_grid, ...
                       dt_window, dt_fft, fftoverlap, ffttaper);
                         
%
[Suu, ~, ~, dof_u] = compute_cpsd(u, u, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
Svv = compute_cpsd(v, v, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
%
Sup = compute_cpsd(u, p, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
Svp = compute_cpsd(v, p, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);
%
Suv = compute_cpsd(u, v, dtime_data, dtime_grid, dt_window, dt_fft, fftoverlap, ffttaper);


%% From pressure, ADCP (u, v)

% Direction where waves come from
% a1 = -real(Sup) ./ (sqrt(Spp.*(Suu + Svv)));
% b1 = -real(Svp) ./ (sqrt(Spp.*(Suu + Svv)));
% Direction where waves go to
a1 = real(Sup) ./ (sqrt(Spp.*(Suu + Svv)));
b1 = real(Svp) ./ (sqrt(Spp.*(Suu + Svv)));
%
a2 = (Suu - Svv)./(Suu + Svv);
b2 = (2.*real(Suv)) ./ (Suu + Svv);


