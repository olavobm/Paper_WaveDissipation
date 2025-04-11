%% Export L2 data files to L3 files, where they all have
% the same time vector and instruments that were swapped
% are merged into a single L3 for that location.

clear
close all


%% List of all sites (same as in the L2 processing)

% All sites
list_sites = ["B01", "B03", "B04", "B05", "B06", "B07", ...
              "B08", "B09", "B10", "B11", "B12", "B13", ...
              "B14", "B15", "B16", "B17", "B18", ...
              "A01", "A02", "A04", "A05", "A06", "A07", "A09", ...
              "C01", "C02", "C05", "C06", "C07", "C08", "C09", ...
              "D01", "D02", ...
              "E01", "E02", "E05", "E07", "E08", "E09", "E11", "E13", ...
              "E03", "E04", "E06", "E12", ...
              "X01", "X03", "X04", "X05", "X06", "X07", ...
              "X08", "X09", "X10", "X11", "X12", "X13", "X14"];


%% Select instruments (matching folder names)

%
list_Instrument_folders = ["ADCPs" "SoloD", "Spotter"];


%% Directory of Level2 and Level3 data

%
dir_L2 = fullfile(paper_directory(), 'data', 'level_2');

%
dir_outputL3 = fullfile(paper_directory(), 'data', 'level_3', 'per_location');


%% Set time limits of the grid

% % %
% % dtime_commonlims = [datetime(2022, 06, 15, 20, 00, 00), ...
% %                     datetime(2022, 07, 21, 05, 00, 00)];
           
%
dtime_grid = get_parameters_procL2();

%
dtime_commonlims = dtime_grid([1, end]);


%%
% ------------------------------------------------------------------
% ------------------ GET L2 DATA AND EXPORT TO L3 ------------------
% ------------------------------------------------------------------


%% Loop over all files, load data, and export to L3

%
for i1 = 1:length(list_sites)

    %
    llookfor_Instrument = true;
    ind_lookfor = 0;
    
    %
    while llookfor_Instrument

        %% Find files ...

        %
        ind_lookfor = ind_lookfor + 1;
        
        %
        dir_L2aux = fullfile(dir_L2, [char(list_Instrument_folders(ind_lookfor)) '_Level2']);

        %
        alldata_files = dir(fullfile(dir_L2aux, '*.mat'));
        
        
        %% ... that match site (can be multiple files per site)
        
        %
        ind_matchID = [];
        
        %
        for i2 = 1:length(alldata_files)
            %
            ind_aux = strfind(alldata_files(i2).name, list_sites(i1));
            %
            if ~isempty(ind_aux)
                %
                ind_matchID = [ind_matchID, i2];
            end
        end

        % If could not find ind_matchID, then go to next instrument type
        if ~isempty(ind_matchID)
            llookfor_Instrument = false;
        end
        
    end

    
    %% Load data

    %
    if length(ind_matchID)==1

        %
        data_aux = load(fullfile(dir_L2aux, alldata_files(ind_matchID).name));
        field_aux = fieldnames(data_aux); 
        %
        data_aux = data_aux.(field_aux{1});

    %
    elseif length(ind_matchID)==2 && (strcmp(list_sites(i1), "B01") || ...
                                      strcmp(list_sites(i1), "E07"))

        %
        data_aux = Spotter_mergeL2(alldata_files(ind_matchID(1)).name, ...
                                   alldata_files(ind_matchID(2)).name);

    %
    else

        error('unexpected number of files.')

    end


    %% Time-rebase data to common time grid
    % (i.e., for timeseries shorter than the time common time grid,
    % pad the edges with NaNs so that timeseries from all sites
    % can be put in the same array)

    %
    dataL3 = roxsi_retimeL2(data_aux, dtime_commonlims);

    
    %% Run this block for finding the time
    % limits that encompasses all the dataset
    %
    % otherwise, comment it out.
    
% %     %
% %     fields_aux = fieldnames(dataL3);
% %     
% %     %
% %     if any(strcmp(fields_aux, 'bottomdepth'))
% %         str_bottom = 'bottomdepth';
% %     elseif any(strcmp(fields_aux, 'bottomdepthfrompres'))
% %         str_bottom = 'bottomdepthfrompres';
% %     else
% %         error('variable not found')
% %     end
% %     
% %     
% %     %
% %     if i1==1
% %         
% %         %
% %         ind_1 = find(~isnan(dataL3.(str_bottom)), 1, 'first');
% %         ind_2 = find(~isnan(dataL3.(str_bottom)), 1, 'last');
% %         %
% %         dtime_pick = dataL3.dtime([ind_1, ind_2]);
% %         
% % 	%
% %     else
% % 
% %         %
% %         ind_1 = find(~isnan(dataL3.(str_bottom)), 1, 'first');
% %         ind_2 = find(~isnan(dataL3.(str_bottom)), 1, 'last');
% %         
% %         %
% %         if dataL3.dtime(ind_1) < dtime_pick(1)
% %             dtime_pick(1) = dataL3.dtime(ind_1);
% %         end
% %         %
% %         if dataL3.dtime(ind_2) > dtime_pick(2)
% %             dtime_pick(2) = dataL3.dtime(ind_2);
% %         end
% %         
% %     end
        

    %% Save to L3 data folder

    %
    save(fullfile(dir_outputL3, ['roxsi_L3_' char(list_sites(i1)) '.mat']), "dataL3")
    

    
end


%% Progress message


disp('-------- DONE CREATING L3 DATA FILES AT EACH SITE --------')
