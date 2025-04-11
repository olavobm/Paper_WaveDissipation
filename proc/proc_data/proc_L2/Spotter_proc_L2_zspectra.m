function [frequency, Szz, indgridlims] = Spotter_proc_L2_zspectra(spotterL1, dtime_grid, fftparams)
%% [frequency, Szz, indgridlims] = SPOTTER_PROC_L2_ZSPECTRA(spotterL1, dtime_grid, fftparams)
%
%   inputs
%       - spotterL1: L1 ROXSI Spotter data structure.
%       - dtime_grid: datetime vector.
%       - fftparams: structure variable with parameters for spectra.
%
%   outputs
%       - frequency: frequency vector
%       - Szz: matrix with sea surface elevation spectra.
%       - indgridlims: 1x2
%
%
% SPOTTER_PROC_L2_ZSPECTRA.m is a high-level function to compute sea
% surface elevation spectra from a L1 Spotter data structure in ROXSI.
% The spectra is computed by a lower-level function that uses Matlab's
% cpsd.
%
%
% See also:
%   compute_cpsd.m



%%

%
[Szz_tmp, frequency, indgridlims] = ...
            compute_cpsd(spotterL1.displacement.z, spotterL1.displacement.z, ...
                         spotterL1.displacement.dtime, dtime_grid, ...
                         fftparams.dt_avgwindow, fftparams.dt_fft, ...
                         fftparams.fftoverlap, fftparams.ffttaper);


%% Remove unreliable estimates a low frequencies.
% (In this block code, keep spectra only between initial
% deployment and final recovery (though with NaNs in
% between if gaps exist))
  
%
Szz = Szz_tmp;

%
lbelowcutoff = (frequency <= fftparams.lowfreq_cutoff);
%
Szz(:, lbelowcutoff) = NaN;


%% Keeping the spectra with for all grid points, with NaNs
% before/after data is available
            
% % %
% % Szz = NaN(length(dtime_grid), length(frequency));              
% %                      
% % %
% % ind_fill = indspeclims(1) : 1 : indspeclims(2);
% % 
% % %
% % Szz(ind_fill, :) = Szz_tmp;
% % 
% % %
% % lbelowcutoff = (frequency <= fftparams.lowfreq_cutoff);
% % 
% % %
% % Szz(:, lbelowcutoff) = NaN;




