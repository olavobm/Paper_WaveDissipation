function adcp_proc_L2(dirdataL1, filename, dirdataL2, dtime_grid, fftparams)
%% ADCP_PROC_L2(dirdataL1, filename, dirdataL2, dtime_grid, fftparams)
%
%   inputs
%       - dirdataL1: directory with L1 ADCP data.
%       - filename: name of one file with L1 ADCP data.
%       - dirdataL2: directory where L2 ADCP data will be saved.
%       - dtime_grid: datetime vector for spectral estimates.
%       - fftparams: structure with parameters for computing spectra.
%
%
% ADCP_PROC_L2.m is a high-level function for doing the L1-to-L2
% data processing of ADCP data from the ROXSI 2022 experiment.
%
% Fields of fftparams structure:
%       (dt_fft, fftoverlap, ffttaper, dt_avgwindow, lowfreq_cutoff)
%
%
% See also:
%   compute_binaverage.m
%   compute_cpsd.m
%   adcp_proc_dirmoments.m
%   wave_meandir_dirspread.m
%   wave_freqtok.m
%   wave_cp.m
%   wave_cg.m



%%

tic

%
str_filename_input = char(filename);
%
ind_match_A_aux = strfind(str_filename_input, '_L1_');
ind_match_B_aux = strfind(str_filename_input, '.mat');
%
ind_get_aux = (ind_match_A_aux + 4) : (ind_match_B_aux-1);

%
disp(['--- Processing ADCP L2 data: ' str_filename_input(ind_get_aux)])


%%

%
freqbounds_to_k = [0.04, fftparams.highfreq_cutoff];


%% Load L1 velocity data

%
disp('Loading L1 data ...')

%
dataL1 = load(fullfile(dirdataL1, filename));
field_aux = fieldnames(dataL1);
dataL1 = dataL1.(field_aux{1});


disp('Done loading data')


%% Create and adcpL2 by copying metadata from dataL1

%
adcpL2.locationID = dataL1.locationID;
adcpL2.SN = dataL1.SN;
adcpL2.instrument = dataL1.instrument;
%
adcpL2.latitude = dataL1.latitude;
adcpL2.longitude = dataL1.longitude;
adcpL2.site = dataL1.site;
adcpL2.X = dataL1.X;
adcpL2.Y = dataL1.Y;

%
adcpL2.rho0 = dataL1.rho0;
adcpL2.g = dataL1.g;


%
adcpL2.zhab_transducer = dataL1.zhab_transducer;
%
adcpL2.zhab_vel = dataL1.zhab_vel;


%% Add FFT parameters to adcpL2

% Parameters for fft/spectra
adcpL2.fftparams = fftparams;


%% Pass the time vector of the spectra to adcpL2

%
adcpL2.dtime = dtime_grid(:);


%% Compute bottom depth (hourly average quantities from L1 data)

%
adcpL2.bottomdepthfrompres = compute_binaverage(dataL1.bottomdepthfrompres, dataL1.dtime, ...
                                                adcpL2.dtime, diff(adcpL2.dtime(1:2)));
%
adcpL2.bottomdepthfrompres = adcpL2.bottomdepthfrompres(:);


%% Compute pressure spectra

%
disp('--- Computing spectra ---')

% Compute spectra
[Spp_tmp, frequency, ...
 indgrid_bounds, dof_aux] = ...
            compute_cpsd(dataL1.pressure, dataL1.pressure, ...
                         dataL1.dtime, dtime_grid, ...
                         fftparams.dt_avgwindow, fftparams.dt_fft, ...
                         fftparams.fftoverlap, fftparams.ffttaper);

%
lbelowcutoff = (frequency <= fftparams.lowfreq_cutoff);
%
Spp_tmp(:, lbelowcutoff) = NaN;

% Trim out high-frequencies
lbelow_highcutoff = (frequency <= fftparams.highfreq_cutoff);

%
frequency = frequency(lbelow_highcutoff);
Spp_tmp = Spp_tmp(:, lbelow_highcutoff);

%
adcpL2.frequency = frequency;
%
adcpL2.Spp = Spp_tmp;

%
adcpL2.dof_pp = dof_aux;

% Trim time for times when data is available
adcpL2.dtime = adcpL2.dtime(indgrid_bounds(1):indgrid_bounds(2));

%
disp('--- Spectra computed ---')


%% Compute directional moments (a1, b1, a2, b2) and
% the 2 sets of mean direction and directional spread

%
disp('--- Computing directional moments ---')

%
[adcpL2.a1, adcpL2.b1, ...
 adcpL2.a2, adcpL2.b2, ~, ~, ~, ~, ...
 dof_u] = adcp_proc_dirmoments(dataL1, adcpL2.dtime, adcpL2.fftparams);

%
adcpL2.dof_uu = dof_u;

%
disp('--- Done with directional moments ---')



%% Compute wave direction and directional spread

%
[adcpL2.meandir1, adcpL2.dirspread1, ...
 adcpL2.meandir2, adcpL2.dirspread2] = wave_meandir_dirspread(adcpL2.a1, adcpL2.b1, ...
                                                              adcpL2.a2, adcpL2.b2);

% 
adcpL2.meandir1 = (180/pi) * adcpL2.meandir1;
adcpL2.meandir2 = (180/pi) * adcpL2.meandir2;
%
adcpL2.dirspread1 = (180/pi) * adcpL2.dirspread1;
adcpL2.dirspread2 = (180/pi) * adcpL2.dirspread2;

% % %  if meandir1 was direction of where waves are coming from, then
% % %  do this to transform towards direction of wave propagation.
% % adcpL2.meandir1 = (180/pi) * wrapToPi(adcpL2.meandir1 + pi);



%% Compute wavenumber

%
lfreqinlims = (adcpL2.frequency >= freqbounds_to_k(1)) & ...
              (adcpL2.frequency <= freqbounds_to_k(2));

%
adcpL2.k = NaN(length(adcpL2.dtime), length(adcpL2.frequency));

%
for i = 1:length(adcpL2.dtime)
    
    %
    adcpL2.k(i, lfreqinlims) = wave_freqtok(adcpL2.frequency(lfreqinlims), ...
                                            adcpL2.bottomdepthfrompres(i));

end


%% Compute transfer function

%
adcpL2.K2 = (cosh(adcpL2.k .* adcpL2.bottomdepthfrompres) ./ ...
             cosh(adcpL2.k * adcpL2.zhab_transducer)).^2;


%% Compute surface elevation spectra

adcpL2.See = adcpL2.Spp .* adcpL2.K2;


%% Compute phase speed and group velocity

adcpL2.cp = wave_cp(adcpL2.k, adcpL2.bottomdepthfrompres);
adcpL2.cg = wave_cg(adcpL2.k, adcpL2.bottomdepthfrompres);


%% Compute energy flux through linear wave theory

%
adcpL2.E = (adcpL2.rho0 * adcpL2.g) * adcpL2.See;
%
adcpL2.Ecg = adcpL2.E .* adcpL2.cg;

    
%% Compute cross-shore and alongshore NET flux with moments
%
% (because Ecg is NOT a NET Flux, Fx and Fy are NOT
% components of the total flux).

%
adcpL2.Fx = adcpL2.Ecg .* adcpL2.a1;
adcpL2.Fy = adcpL2.Ecg .* adcpL2.b1;


%%
% --------------------------------------------------
% ----------- SAVE THE L2 DATA STRUCTURE -----------
% --------------------------------------------------


%% Add README

%
dtime_proc = datetime('now', 'TimeZone', 'Local');
%
time_dataproc_char = datestr(dtime_proc, 'yyyy/mm/dd HH:MM:SS');

% Add a README
adcpL2.README = ['Level 2 ADCP data from the ROXSI 2022 large-scale ' ...
                'array. The data is from an ADCP with serial number ' ...
                'SN and deployed at instrument site locationID. ' ...
                'Data processed by script ' mfilename() '.m on ' ...
                time_dataproc_char '. Net flux components Fx and Fy ' ...
                'are in the local XY coordinate system. meandir1 is ' ...
                'direction where waves are propagating to, in ' ...
                'counterclockwise degrees from +X.'];


%% Save L2 data file

%
disp('----- Saving level 2 data -----')
str_filename = ['roxsi_adcp_L2_' char(adcpL2.locationID) '_' char(adcpL2.SN) '.mat'];
%
save(fullfile(dirdataL2, str_filename), 'adcpL2', '-v7.3')
%
toc
