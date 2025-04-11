function Spotter_proc_b_L2(dirpaper, dir_spotterL2, ID_SN, strdepth, filedepth)
%% SPOTTER_PROC_B_L2(dirpaper, dir_spotterL2, ID_SN, strdepth, filedepth)
%
%   inputs
%       - dirpaper: parent directory of the paper.
%       - dir_spotterL2Ç directory with L2 Spotter data.
%       - ID_SN: characters with mooring ID and Spotter serial number.
%       - strdepth: 'bathy' or 'pressure'.
%       - filedepth (optional): file with bathymetry data.
%
%
% SPOTTER_PROC_b_L2.m is a high-level function that incorporates bathymetry
% and tide data into the L2 Spotter data. The bottom depth is then used
% to calculate statistics that depend on bottom depth (e.g., group
% velocity, energy flux, etc).
%
%
% See also:
%   run_proc_Spotter_L2.m
%   Spotter_proc_a_L2.m
%   wave_freqtok.m
%   wave_cp.m
%   wave_cg.m


%% Define files to get the data from
%
% The files from the effective depth are used, but using the smoothed
% depth from bathy_longscale_*.mat gives very similar results. Note
% that the effective depth correction of the spectra is only used
% for pressure sensor data, not for Spotter data.

%
filebathy.Asilomar = fullfile(dirpaper, 'data', 'bathymetry', 'bathy_dhbathy_Asilomar.mat');
filebathy.ChinaRock = fullfile(dirpaper, 'data', 'bathymetry', 'bathy_dhbathy_ChinaRock.mat');

%
if exist('filedepth', 'var')
    ldefaultbathy = false;
else
    ldefaultbathy = true;
end

%
filetide = fullfile(dirpaper, 'data', 'level_1', 'NOAA_data', 'noaa_tides.mat');


%% Frequency bound to compute wavenumber

freqbounds_to_k = [0.04, 1.1];


%% Print progress message on the screen

%
disp(['--- Processing L1-to-L2 (part B) Spotter data: ' ID_SN ' ---'])


%% Load L2 data

%
filename_L2data = ['roxsi_spotter_L2_' ID_SN '.mat'];

%
spotterL2 = load(fullfile(dir_spotterL2, filename_L2data));
spotterL2 = spotterL2.spotterL2;


%% Get bottom depth
% 
% could be a co-located pressure sensor
% or requires bathymetry + tides
% or requires some large-scale bathymetry (Spotter P01 offshore) + tide
    
%
if strcmp(strdepth, 'bathy')

    %
    if ldefaultbathy
        bathydata = load(filebathy.(spotterL2.site));
    else
        bathydata = load(filedepth);
    end
    
    %
    field_aux = fieldnames(bathydata);
    %
    if length(field_aux)==1
        bathydata = bathydata.(field_aux{1});
    end
    
    %
    tidedata = load(filetide);
    field_aux = fieldnames(tidedata);
    %
    if length(field_aux)==1
        tidedata = tidedata.(field_aux{1});
    end
    
    
    % ------------------------------
    % Removing NaNs
    
    %
    [xgrid, ygrid] = meshgrid(bathydata.x, bathydata.y);
    %
    xgrid = xgrid(:);
    ygrid = ygrid(:);
    z_msl = bathydata.zmsl(:);
    %
    lokbathy = ~isnan(z_msl);
    
    % 
    interp_zmsl = scatteredInterpolant(xgrid(lokbathy), ygrid(lokbathy), ...
                                       z_msl(lokbathy));
	%
    spotterL2.location.zmsl = interp_zmsl(spotterL2.location.X, ...
                                           spotterL2.location.Y);
    
	% ------------------------------
	% Add tides to the MSL bathymetry    
    
    %
    spotterL2.location.z_tide = interp1(tidedata.dtime, tidedata.z, spotterL2.dtime);
    
    %
    spotterL2.bottomdepth = -spotterL2.location.zmsl + ...
                             spotterL2.location.z_tide;
                         
end


%% Compute wavenumber

%
lfreqinlims = (spotterL2.frequency >= freqbounds_to_k(1)) & ...
              (spotterL2.frequency <= freqbounds_to_k(2));

%
spotterL2.k = NaN(length(spotterL2.dtime), length(spotterL2.frequency));

%
for i = 1:length(spotterL2.dtime)
    %
    spotterL2.k(i, lfreqinlims) = wave_freqtok(spotterL2.frequency(lfreqinlims), ...
                                               spotterL2.bottomdepth(i));
    
end


%% Compute phase speed and group velocity

spotterL2.cp = wave_cp(spotterL2.k, spotterL2.bottomdepth);
spotterL2.cg = wave_cg(spotterL2.k, spotterL2.bottomdepth);


%% Compute energy flux through linear wave theory
    
%
spotterL2.rho0 = 1025;
spotterL2.g = 9.8;

%
spotterL2.E = (spotterL2.rho0 * spotterL2.g) * spotterL2.See;
%
spotterL2.Ecg = spotterL2.E .* spotterL2.cg;

    
%% Compute cross-shore and alongshore NET fluxes
% with first directional moments. These are net fluxes
% of where waves are going to.
%
% (because Ecg is NOT a NET Flux, Fx and Fy are NOT
% components of the total flux).

%
spotterL2.Fx = spotterL2.Ecg .* spotterL2.a1;
spotterL2.Fy = spotterL2.Ecg .* spotterL2.b1;


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
spotterL2.README = ['Level 2 Spotter data from ROXSI 2022 large-scale ' ...
                    'array. The data is from Spotter with serial number ' ...
                    'SN and deployed at instrument site locationID. ' ...
                    'Data processed by script ' mfilename() '.m on ' ...
                    time_dataproc_char '. Horizontal displacements ' ...
                    '(x, y) are in the local XY coordinate system.'];


%% Save/overwrite L2 data file

%
str_filename = ['roxsi_spotter_L2_' char(spotterL2.locationID) '_' char(spotterL2.SN) ''];
%
save(fullfile(dir_spotterL2, [str_filename '.mat']), 'spotterL2', '-v7.3')
%
disp('--- Level 2 data saved ---')
disp(['--- DONE with L1-to-L2 (part B) Spotter data processing: ' ID_SN ' ---'])




