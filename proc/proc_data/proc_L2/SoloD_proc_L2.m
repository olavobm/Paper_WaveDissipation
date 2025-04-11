function SoloD_proc_L2(dirdataL1, dirdataL2, ID_SN, dtime_grid, fftparams)
%% SOLOD_PROC_L2(dirdataL1, dirdataL2, ID_SN, dtime_grid, fftparams)
%
%   inputs
%       - dirdataL1: directory with L1 data.
%       - dirdataL2: directory with L2 data.
%       - ID_SN: characters with mooring ID and SoloD serial number.
%       - dtime_grid: time grid for computing spectra.
%       - fftparams: structure where fields have parameters for
%                    computing spectra.
%
%
% SOLOD_PROC_L2.m is a high-level function for data processing where
% spectra are computed from SoloD data. Results (i.e. L2 data) are saved
% in the directory dirdataL2.
%
% Fields of fftparams structure:
%       (dt_fft, fftoverlap, ffttaper, dt_avgwindow, lowfreq_cutoff)
%
%
% See also:
%   run_proc_SoloD_L2.m
%   get_parameters_procL2.m
%   compute_cpsd.m
%   compute_binaverage.m
%   wave_freqtok.m


%% Display progress message

tic

%
disp(['Processing SoloD: ' ID_SN])


%% Load L1 data

%
filename_L1data = ['roxsi_soloD_L1_' ID_SN(1:3) '_' ID_SN(5:end) '.mat'];

%
solodL1 = load(fullfile(dirdataL1, filename_L1data));
solodL1 = solodL1.solodL1;


%% Create and soloDL2 by copying metadata from spotterL1

%
solodL2.locationID = solodL1.locationID;
solodL2.SN = solodL1.SN;
%
solodL2.instrument = "SoloD";
%
solodL2.latitude = solodL1.latitude;
solodL2.longitude = solodL1.longitude;
solodL2.site = solodL1.site;
solodL2.X = solodL1.X;
solodL2.Y = solodL1.Y;
%
solodL2.zhab = solodL1.zhab;


%% Normalize pressure to have it in units of meter

%
solodL2.rho0 = 1025;
solodL2.g = 9.8;

% Go from dbar to m
solodL1.pressure = (1e4*solodL1.pressure)./(solodL2.rho0*solodL2.g);


%% Add FFT parameters to soloDL2

% Parameters for fft/spectra
solodL2.fftparams = fftparams;


%% Frequency bounds to compute wavenumber

%
freqbounds_to_k = [0.04, 1];


%% Pass the time vector of the spectra to soloDL2

solodL2.dtime = dtime_grid(:);


%% Compute pressure spectra

%
[Spp_tmp, frequency, indgrid_bounds, dof] = ...
            compute_cpsd(solodL1.pressure, solodL1.pressure, ...
                         solodL1.dtime, dtime_grid, ...
                         fftparams.dt_avgwindow, fftparams.dt_fft, ...
                         fftparams.fftoverlap, fftparams.ffttaper);

%
solodL2.dtime = solodL2.dtime(indgrid_bounds(1):indgrid_bounds(2));
%
solodL2.frequency = frequency;


% Remove unreliable estimates a low frequencies.
% (In this block code, keep spectra only between initial
% deployment and final recovery (though with NaNs in
% between if gaps exist))
%
Spp = Spp_tmp;
%
lbelowcutoff = (frequency <= fftparams.lowfreq_cutoff);
%
Spp(:, lbelowcutoff) = NaN;
%
solodL2.Spp = Spp;


%% Compute averaged bottom depth from pressure

%
[p_mean, indgridlims] = compute_binaverage(solodL1.pressure, solodL1.dtime, ...
                                           solodL2.dtime, fftparams.dt_avgwindow);

%
solodL2.dtime = solodL2.dtime(indgridlims(1):indgridlims(2));
%
solodL2.bottomdepthfrompres = p_mean(:);
%
solodL2.bottomdepthfrompres = solodL2.bottomdepthfrompres + solodL2.zhab;


%% Compute wavenumber

%
lfreqinlims = (solodL2.frequency >= freqbounds_to_k(1)) & ...
              (solodL2.frequency <= freqbounds_to_k(2));

%
solodL2.k = NaN(length(solodL2.dtime), length(solodL2.frequency));

%
for i = 1:length(solodL2.dtime)
    
    %
    solodL2.k(i, lfreqinlims) = wave_freqtok(solodL2.frequency(lfreqinlims), ...
                                             solodL2.bottomdepthfrompres(i));
    
end


%% Compute transfer function

%
h_array = repmat(solodL2.bottomdepthfrompres, 1, size(solodL2.k, 2));

%
solodL2.K2 = (cosh(solodL2.k.*h_array) ./ ...
              cosh(solodL2.k * solodL2.zhab)).^2;


%% Compute surface elevation spectra

solodL2.See = solodL2.Spp .* solodL2.K2;


%% Compute phase speed and group velocity

solodL2.cp = wave_cp(solodL2.k, solodL2.bottomdepthfrompres);
solodL2.cg = wave_cg(solodL2.k, solodL2.bottomdepthfrompres);


%% Compute energy flux through linear wave theory
    
%
solodL2.E = (solodL2.rho0 * solodL2.g) * solodL2.See;
%
solodL2.Ecg = solodL2.E .* solodL2.cg;


%%
% --------------------------------------------------
% ----------- SAVE THE L2 DATA STRUCTURE -----------
% --------------------------------------------------

%% Save L2 data file

%
disp('----- Saving level 2 data -----')
str_filename = ['roxsi_soloD_L2_' char(solodL2.locationID) '_' char(solodL2.SN) ''];
%
save(fullfile(dirdataL2, [str_filename '.mat']), 'solodL2', '-v7.3')
%
disp(['----- DONE with L1 to L2 SoloD data processing: ' ID_SN ' -----'])

%
disp(' ')
disp('----- Run processing time was:')
%
toc
%
disp(' ')







