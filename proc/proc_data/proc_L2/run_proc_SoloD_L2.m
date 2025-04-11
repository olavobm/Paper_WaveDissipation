%% High-level script that processes SoloD data (from Level 1 to
% Level 2) from the large-scale array in the ROXSI's 2022 experiment.

clear
close all


%%

%
dir_dataL1 = fullfile(paper_directory(), 'data', 'level_1', 'SoloD_Level1');
%
dir_dataL2 = fullfile(paper_directory(), 'data', 'level_2', 'SoloD_Level2');


%% Choose SoloDs to process data

%
list_SoloD = {'A02_077589', 'A04_077590', 'A05_077817', ...
              'A06_077587', 'A07_077848', 'A09_077591', ...
...%
              'B04_124020', 'B06_124022', 'B07_124073', ...
              'B09_077845', 'B12_077813', 'B14_077808', ...
              'B16_077586', 'B18_077272', ...
...%
              'C02_077594', 'C03_124039', 'C05_077595', ...
              'C06_077809', 'C07_077819', 'C08_077804', 'C09_077274', ...
...%
              'D01_124015', 'D02_124017', ...
...%
              'X07_077810', 'X08_077818', 'X09_077814', ...
              'X10_077807', 'X12_077815', 'X14_077816'};

          
%% Get common parameters for L2 data processing

%
[dtime_grid, fftparams] = get_parameters_procL2();


%% Print progress message on the screen

%
totalRunTime = tic;

%
disp(' '), disp(' ')
disp('----- Processing SoloD: L1 to L2 -----')
disp(['Started running the code at: ' datestr(now())])


%% Loop over SoloDs and do data processing for each


% Looping over Spotters
for i = 1:length(list_SoloD)

    
    %% Run L2 data processing that only requires Spotter data
    
    %
    SoloD_proc_L2(dir_dataL1, dir_dataL2, list_SoloD{i}, dtime_grid, fftparams)
    
    
    %% Print progress message
    
    disp(['Done with L2 processing of SoloD ' list_SoloD{i}])
    
end



%% Final progress message

%
disp('--- DONE with L2 processing for all SoloDs ---')
disp('--- Run processing time was:')
%
toc(totalRunTime)


