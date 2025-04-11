%% This is just the same as compute_Suu_SAb.m and compute_bulkstatistics.m,
% but doing for the file that has the L3 data with local depth hp.

clear
close all


%% Load L3 data

%
dir_dataL3 = fullfile(paper_directory(), 'data', 'level_3');

%
dataL3 = load(fullfile(dir_dataL3, 'roxsi_L3_allsites_hp.mat'));
dataL3 = dataL3.dataL3;

%
dataL4 = dataL3;


%% Compute spectra of U bottom and A bottom

%
freq_array = reshape(dataL4.frequency, [1, length(dataL4.frequency), 1]);
%
h_array = reshape(dataL4.bottomdepth, [length(dataL4.dtime), 1, length(dataL4.locationID)]);

%
dataL4.Suu = dataL4.See .* (2*pi*freq_array ./ sinh(dataL4.k .* h_array)).^2;
%
dataL4.SAb = dataL4.See .* (1 ./ sinh(dataL4.k .* h_array)).^2;


%% Define frequency bands for integrating the spectra

%
list_freqbands = "seaswell";
freq_lims = [0.05, 0.2];


%% Compute sea-swell bulk properties

%
linfreq = (dataL4.frequency >= freq_lims(1)) & ...
          (dataL4.frequency <= freq_lims(2));
      
%
dataL4.list_freqbands = list_freqbands;

%
dataL4.(list_freqbands).freqlims = freq_lims;


%% Compute Hs and mean frequency

%
m0 = trapz(dataL4.frequency(linfreq), dataL4.See(:, linfreq, :), 2);
m0 = squeeze(m0);

%
m1 = trapz(dataL4.frequency(linfreq), dataL4.See(:, linfreq, :) .* ...
                                      dataL4.frequency(linfreq).', 2);
m1 = squeeze(m1);

%
dataL4.(list_freqbands).Hs = 4*sqrt(m0);

%
dataL4.(list_freqbands).meanfrequency = m1./m0;
    
    
%% Compute Urms and Ab

%
dataL4.(list_freqbands).Urms = ...
                    sqrt(trapz(dataL4.frequency(linfreq), ...
                               dataL4.Suu(:, linfreq, :), 2));
dataL4.(list_freqbands).Urms = squeeze(dataL4.(list_freqbands).Urms);

%
dataL4.(list_freqbands).Ab = ...
                    sqrt(2)*sqrt(trapz(dataL4.frequency(linfreq), ...
                                       dataL4.SAb(:, linfreq, :), 2));
dataL4.(list_freqbands).Ab = squeeze(dataL4.(list_freqbands).Ab);


%% Compute directional moments

list_moments = ["a1", "b1", "a2", "b2"];

%
for i = 1:length(list_moments)
    
    %
    var_aux = trapz(dataL4.frequency(linfreq), ...
                    dataL4.See(:, linfreq, :) .* ...
                    dataL4.(list_moments(i))(:, linfreq, :), 2);
    var_aux = squeeze(var_aux);
    
	%
    dataL4.(list_freqbands).(list_moments(i)) = var_aux ./ m0;
    
end


%% Compute mean direction and directional spread

%
[dataL4.(list_freqbands).meandir1, ...
 dataL4.(list_freqbands).dirspread1, ...
 dataL4.(list_freqbands).meandir2, ...
 dataL4.(list_freqbands).dirspread2] = ...
                    wave_meandir_dirspread(dataL4.(list_freqbands).a1, ...
                                           dataL4.(list_freqbands).b1, ...
                                           dataL4.(list_freqbands).a2, ...
                                           dataL4.(list_freqbands).b2);

% 
% % sigL2.meandir1 = (180/pi) * wrapToPi(sigL2.meandir1 + pi);  % if this is direction of where waves coming from
dataL4.(list_freqbands).meandir1 = (180/pi) * dataL4.(list_freqbands).meandir1;   % if this is direction of where waves are going to
dataL4.(list_freqbands).meandir2 = (180/pi) * dataL4.(list_freqbands).meandir2;
%
dataL4.(list_freqbands).dirspread1 = (180/pi) * dataL4.(list_freqbands).dirspread1;
dataL4.(list_freqbands).dirspread2 = (180/pi) * dataL4.(list_freqbands).dirspread2;


%% Compute integrated Ecg

%
dataL4.(list_freqbands).Ecg = squeeze(trapz(dataL4.frequency(linfreq), ...
                                            dataL4.Ecg(:, linfreq, :), 2));
    

%% Compute integrated Fx and Fy

%
dataL4.(list_freqbands).Fx = squeeze(trapz(dataL4.frequency(linfreq), ...
                                           dataL4.Fx(:, linfreq, :), 2));
%
dataL4.(list_freqbands).Fy = squeeze(trapz(dataL4.frequency(linfreq), ...
                                           dataL4.Fy(:, linfreq, :), 2));
                                

%% Save L4 data

%
dir_output = fullfile(paper_directory(), 'data', 'level_4');

%
save(fullfile(dir_output, 'roxsi_datahp_L4.mat'), 'dataL4')
