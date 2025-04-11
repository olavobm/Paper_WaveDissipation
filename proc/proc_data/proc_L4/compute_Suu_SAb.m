%% Compute velocity and horizontal excursion spectra.

clear
close all


%% Load L3 data

%
dir_dataL3 = fullfile(paper_directory(), 'data', 'level_3');

%
dataL3 = load(fullfile(dir_dataL3, 'roxsi_L3_allsites.mat'));
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


%% Save (overwrite) L3 data

%
dir_output = fullfile(paper_directory(), 'data', 'level_4');

%
save(fullfile(dir_output, 'roxsi_data_L4.mat'), 'dataL4')
