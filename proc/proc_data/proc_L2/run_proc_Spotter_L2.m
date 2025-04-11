%% High-level script that processes Spotter data (from Level 1 to
% Level 2) from the large-scale array in the ROXSI's 2022 experiment.

clear
close all


%%

dir_dataL1 = fullfile(paper_directory(), 'data', 'level_1', 'Spotter_Level1');
dir_dataL2 = fullfile(paper_directory(), 'data', 'level_2', 'Spotter_Level2');



%% Choose Spotters to process data

%
list_Spotter = {'B01_1158', 'B01_1150', 'B03_1152', 'B05_1153', ...
                'E01_1851', 'E02_1859', 'E05_1853', ...
                'E07_1855', 'E07_1857', ...
                'E08_1852', ...
                'E09_1856', ...
                'E11_1860', 'E13_1849', ...
                'X01_1151', 'X03_1157', 'X04_1155'};


%% Get common time-grid and parameters for spectra

%
[~, fftparams] = get_parameters_procL2();

% For the time grid, actually use this other grid that excludes the edges
dtime_lims = [datetime(2022, 06, 15, 20, 00, 00), ...   
              datetime(2022, 07, 20, 05, 00, 00)];
%
dtime_lims.TimeZone = 'America/Los_Angeles';

%
dt_grid = hours(1);
dtime_grid = dtime_lims(1) : dt_grid : dtime_lims(2);


%% Display progress message on the screen:

%
totalRunTime = tic;

%
disp(' '), disp(' ')
disp('--- Processing Spotter data: from L1 to L2 ---')
disp('List of Spotters being processed:')
%
for i = 1:length(list_Spotter)
    disp([num2str(i) ' - ' list_Spotter{i}])
end

%
disp(' '), disp(' ')
% 
disp(['Started running the code at: ' datestr(now())])
disp(' '), disp(' ')


%% Loop over Spotters and do data processing for each

% Looping over Spotters
for i = 1:length(list_Spotter)

    
    %% Run L2 data processing that only requires Spotter data
    
    %
    Spotter_proc_a_L2(dir_dataL1, dir_dataL2, list_Spotter{i}, dtime_grid, fftparams)
    

    %% Run L2 data processing that requires additional data (e.g. water depth)
    
    %
    Spotter_proc_b_L2(paper_directory(), dir_dataL2, list_Spotter{i}, 'bathy')
    
    
    %% Print progress message
    
    %
    disp(['Done with L2 processing of Spotter ' list_Spotter{i}])
    toc(totalRunTime)
    %
    disp(' ')
    
end


%%

%
disp(' ')
disp('#### Done with L2 data processing for all Spotters #####')

%
disp(' ')
disp('*** The total time to run the data processing was:')
%
toc(totalRunTime)











