function Spotter_proc_a_L2(dir_dataL1, dir_dataL2, ID_SN, dtime_grid, fftparams)
%% SPOTTER_PROC_A_L2(dir_dataL1, dir_dataL2, ID_SN, dtime_grid, fftparams)
%
%   inputs
%       - dirdataL1: directory with L1 data.
%       - dirdataL2: directory with L2 data.
%       - ID_SN: characters with mooring ID and Spotter serial number.
%       - dtime_grid: time grid of the spectr.
%       - fftparams: structure with parameters for computing spectra.
%
%
% SPOTTER_PROC_a_L2.m is a high-level function for processing Spotter data
% with spectral calculations. All calculations are based only on the
% measurements from the Spotter wave buoy (as opposed to
% Spotter_proc_b_L2.m).
%
% Fields of fftparams structure:
%       (dt_fft, fftoverlap, ffttaper, dt_avgwindow, lowfreq_cutoff)
%
%
% See also:
%   run_proc_Spotter_L2.m
%   get_parameters_procL2.m
%   Spotter_proc_b_L2.m
%   Spotter_proc_L2_zspectra.m
%   Spotter_proc_dirmoments.m
%   wave_meandir_dirspread.m


%% Print progress message on the screen

%
disp(['--- Processing L1-to-L2 (part A) Spotter data: ' ID_SN ' ---'])


%% Get L1 data file name (spotter or smartmooring_

%
filename_L1data = ['roxsi_spotter_L1_' ID_SN '.mat'];

%
info_file = dir(fullfile(dir_dataL1, filename_L1data));

%
if isempty(info_file)
    filename_L1data = ['roxsi_smartmooring_L1_' ID_SN '.mat'];
end


%% Load L1 data
%
spotterL1 = load(fullfile(dir_dataL1, filename_L1data));
spotterL1 = spotterL1.spotterL1;


%% Create and spotterL2 by copying metadata from spotterL1

%
spotterL2.locationID = spotterL1.locationID;
spotterL2.SN = spotterL1.SN;
spotterL2.instrument = spotterL1.instrument;
%
spotterL2.latitude = spotterL1.latitude;
spotterL2.longitude = spotterL1.longitude;
spotterL2.site = spotterL1.site;
spotterL2.X = spotterL1.X;
spotterL2.Y = spotterL1.Y;


%% Add averaged measured locations

%
[indsub, reshapeNdims, indgridlims] = reshapeintoWindows(spotterL1.location.dtime, dtime_grid);

%
[x_all, y_all] = ROXSI_lltoxy(spotterL1.location.latitude, ...
                              spotterL1.location.longitude, ...
                              spotterL2.site);
 
%
x_avg = mean(reshape(x_all(indsub), reshapeNdims), 1);
y_avg = mean(reshape(y_all(indsub), reshapeNdims), 1);
%
x_avg = x_avg(:);
y_avg = y_avg(:);

%
[spotterL2.location.latitude, ...
 spotterL2.location.longitude] = ROXSI_xytoll(spotterL2.site, x_avg, y_avg);

%
spotterL2.location.X = x_avg;
spotterL2.location.Y = y_avg;
 

%% Add FFT parameters to spotterL2

%
spotterL2.fftparams = fftparams;


%% Pass the time vector of the spectra to spotterL2

spotterL2.dtime = dtime_grid(:);


%% Compute elevation spectra

%
[spotterL2.frequency, ...
 spotterL2.See, ...
 indgrid_bounds] = Spotter_proc_L2_zspectra(spotterL1, spotterL2.dtime, spotterL2.fftparams);

%
spotterL2.dtime = spotterL2.dtime(indgrid_bounds(1):indgrid_bounds(2));

% this indgrid_bounds is not necessary after the full deployment is
% done and all the trimming times for all Spotters are properly defined


%% Compute directional moments (a1, b1, a2, b2) and
% the 2 sets of mean direction and directional spread

%
[spotterL2.a1, spotterL2.b1, ...
 spotterL2.a2, spotterL2.b2] = Spotter_proc_dirmoments(spotterL1, ...
                                                       spotterL2.dtime, ...
                                                       spotterL2.fftparams);


%% Compute wave direction and directional spread

%
[spotterL2.meandir1, spotterL2.dirspread1, ...
 spotterL2.meandir2, spotterL2.dirspread2] = wave_meandir_dirspread(spotterL2.a1, spotterL2.b1, ...
                                                                    spotterL2.a2, spotterL2.b2);

% 
% % spotterL2.meandir1 = (180/pi) * wrapToPi(spotterL2.meandir1 + pi);  % direction that waves propagate to
spotterL2.meandir1 = (180/pi) * spotterL2.meandir1;
spotterL2.meandir2 = (180/pi) * spotterL2.meandir2;
%
spotterL2.dirspread1 = (180/pi) * spotterL2.dirspread1;
spotterL2.dirspread2 = (180/pi) * spotterL2.dirspread2;


%%
% --------------------------------------------------
% ----------- SAVE THE L2 DATA STRUCTURE -----------
% --------------------------------------------------


%% Save L2 data file

%
str_filename = ['roxsi_spotter_L2_' char(spotterL2.locationID) '_' char(spotterL2.SN) ''];
%
save(fullfile(dir_dataL2, [str_filename '.mat']), 'spotterL2', '-v7.3')
%
disp(['--- DONE with L1-to-L2 (part A) Spotter data processing: ' ID_SN ' ---'])




