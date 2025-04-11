%% High-level script that does L1-to-L2 data processing from
% Aquadopps and Signature1000 from the large-scale array in the
% ROXSI 2022 experiment.

clear
close all


%% Set directories

%
dir_dataL1 = fullfile(paper_directory(), 'data', 'level_1', 'ADCPs_Level1');

%
dir_dataoutput = fullfile(paper_directory(), 'data', 'level_2', 'ADCPs_Level2');


%% List of ADCPs to process

%
list_procADCP.locID = ["B08"; "B11"; "E03"; "E04"; "E06"; "E12"; "X06"; "X13"; ...
                       "B10"; "B13"; "B15"; "B17"; "A01"; "C01"; "X05"; "X11"];
                  
list_procADCP.SN = ["13288"; "12280"; "13300"; "13172"; "9736"; "11150"; "13290"; "9945"; ...
                    "103045"; "103046"; "103056"; "101923"; "103043"; "102128"; "100231"; "101941"];
                
list_procADCP.instrument = ["Aquadopp"; "Aquadopp"; "Aquadopp"; "Aquadopp"; ...
                            "Aquadopp"; "Aquadopp"; "Aquadopp"; "Aquadopp"; ...
                            "Signature"; "Signature"; "Signature"; "Signature"; ...
                            "Signature"; "Signature"; "Signature"; "Signature"];
                        
%
Ninstr = length(list_procADCP.locID);


%% Get filenames corresponding to the instruments specified above

%
all_matfiles = dir(fullfile(dir_dataL1, '*.mat'));

%
list_allmatfiles = strings(length(all_matfiles), 1);
%
for i = 1:length(all_matfiles)
    %
    list_allmatfiles(i) = convertCharsToStrings(all_matfiles(i).name);
end


% Get file names corresponding to list_siteID
ind_file = NaN(Ninstr, 1);

%
for i = 1:Ninstr
    
    % Get filename matching with SN
    file_aux = dir(fullfile(dir_dataL1, ['*_' char(list_procADCP.SN(i)) '.mat']));
    
    %
    if length(file_aux)~=1
        error('L1 data from ADCP not found as expected.')
    end
    
    %
    file_name_aux = convertCharsToStrings(file_aux.name);
    
    %
    ind_file(i) = find(strcmp(list_allmatfiles, file_name_aux));
    
end

                
%% Get common time-grid and parameters for spectra

%
[dtime_grid, fftparams] = get_parameters_procL2();


%% Display message on the screen:

%
disp(' '), disp(' ')
disp('--- Processing ADCPs data: from L1 to L2 ---')
disp('List of ADCPs being processed:')
%
for i = 1:Ninstr
    disp([num2str(i) ' - ' char(list_procADCP.instrument(i)) ' at ' ...
                           char(list_procADCP.locID(i)) ' SN = ' char(list_procADCP.SN(i))])
end

%
disp(' '), disp(' ')
% 

%
disp(['Started running the code at: ' datestr(now())])
               
%
totalRunTime = tic;


%% Loop over ADCPs and do data processing for each instrument

%
for i = 1:Ninstr

    %% Run L2 data processing for Signature1000 data    

    %
    disp(' '), disp(' ')
    
    %
    adcp_proc_L2(dir_dataL1, list_allmatfiles(ind_file(i)), dir_dataoutput, ...
                 dtime_grid, fftparams)

    

    %% Print progress message
    
    disp(['Done with L2 processing for ADCP SN-' char(list_procADCP.SN(i)) ' at ' ...
                                                 char(list_procADCP.locID(i))])

    %
    toc(totalRunTime)
    
end


%% Print progress message

%
disp('### DONE WITH L2 DATA PROCESSING FOR ALL ADCPs ###')

%
toc(totalRunTime)


