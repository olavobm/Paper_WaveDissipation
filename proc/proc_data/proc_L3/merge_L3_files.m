%% Script that merges L3 files into one file


clear
close all


%%

freq_lims = [0, 1];


%% All sensors

% All sites
list_locations = ["B01", "B03", "B04", "B05", "B06", "B07", ...
                  "B08", "B09", "B10", "B11", "B12", "B13", ...
                  "B14", "B15", "B16", "B17", "B18", ...
                  "A01", "A02", "A04", "A05", "A06", "A07", "A09", ...
                  "C01", "C02", "C05", "C06", "C07", "C08", "C09", ...
                  "D01", "D02", ...
                  "E01", "E02", "E05", "E07", "E08", "E09", "E11", "E13", ...
                  "E03", "E04", "E06", "E12", ...
                  "X01", "X03", "X04", "X05", "X06", "X07", ...
                  "X08", "X09", "X10", "X11", "X12", "X13", "X14"];


%%

file_name_ID = "allsites";

%%

%
Nlocs = length(list_locations);


%% Set data directory

%
dir_L3 = fullfile(paper_directory(), 'data', 'level_3');
dir_L3_perloc = fullfile(dir_L3, 'per_location');


%%
% -------------------------------------------
% -------------------------------------------
% -------------------------------------------

%%

%
data_aux = load(fullfile(dir_L3_perloc, ['roxsi_L3_' char(list_locations(1)) '.mat']));
data_aux = data_aux.dataL3;
%
Ntime = length(data_aux.dtime);

%
dtime = data_aux.dtime;
freq_vec = data_aux.frequency;
df = diff(data_aux.frequency(1:2));

%
clear data_aux


%%

%
vec_strings = strings(1, Nlocs);
%
vec_number = NaN(Ntime, Nlocs);


%%

%
dataL3.locationID = vec_strings;
dataL3.SN = vec_strings;
dataL3.instrument = vec_strings;
dataL3.site = vec_strings;

%
dataL3.latitude = vec_number;
dataL3.longitude = vec_number;
dataL3.X = vec_number;
dataL3.Y = vec_number;

%
dataL3.fftparameters = [];


%%

%
dataL3.dtime = dtime(:);

%
dataL3.frequency = freq_lims(1) : df : freq_lims(2);
dataL3.frequency = dataL3.frequency(:);

%
dataL3.bottomdepth = vec_number;


%%

%
array_number = NaN(Ntime, length(dataL3.frequency), Nlocs);

%
dataL3.zhab = NaN(Nlocs, 1);
dataL3.zhab_wavedir = NaN(Nlocs, 1);

%
dataL3.Spp = array_number;

%
list_vars_aux = ["See", ...
                 "a1", "b1", "a2", "b2", ...
                 "meandir1", "dirspread1", "meandir2", "dirspread2", ...
                 "k", "cp", "cg", "E", "Ecg", ...
                 "Fx",  "Fy"];
         
%
for i = 1:length(list_vars_aux)
    dataL3.(list_vars_aux(i)) = array_number;
end


%%
% -------------------------------------------
% --- PUT VARIABLES IN THE SAME STRUCTURE ---
% -------------------------------------------


%%

%
disp('---------- Starting to merge L3 files ----------')
%
disp('Merging L3 data from:')
for i = 1:Nlocs
    disp([num2str(i) ' - ' char(list_locations(i))])
end


tic

%%

list_vars_dir = ["a1", "b1", "a2", "b2", ...
                 "meandir1", "dirspread1", "meandir2", "dirspread2", ...
                 "Fx", "Fy"];


%% Loop over files and put all data together

%
for i = 1:Nlocs


    %% Load data structuredataL3
    
    data_aux = load(fullfile(dir_L3_perloc, ['roxsi_L3_' char(list_locations(i)) '.mat']));
    data_aux = data_aux.dataL3;
    
    
    %%
    
    dataL3.locationID(i) = data_aux.locationID;
    dataL3.SN(i) = data_aux.SN;
    dataL3.instrument(i) = data_aux.instrument;
    dataL3.site(i) = data_aux.site;
    
    
    %%
    
    if i==1
        %
        dataL3.fftparameters = data_aux.fftparams;
        %
        dataL3.fftparameters = rmfield(dataL3.fftparameters, ...
                                       {'lowfreq_cutoff', 'highfreq_cutoff'});
    end
    
    
    %%
    
    %
    if strcmp(dataL3.instrument(i), "Spotter")
        %
        dataL3.latitude(:, i) = data_aux.location.latitude;
        dataL3.longitude(:, i) = data_aux.location.longitude;
    else
        %
        dataL3.latitude(:, i) = data_aux.latitude;
        dataL3.longitude(:, i) = data_aux.longitude;
    end
    
    %
    [dataL3.X(:, i), ...
     dataL3.Y(:, i)] = ROXSI_lltoxy(dataL3.latitude(:, i), ...
                                    dataL3.longitude(:, i), ...
                                    dataL3.site(i));
    
    
    %%
    
    %
    l_getfreq = (data_aux.frequency >= dataL3.frequency(1)) & ...
                (data_aux.frequency <= dataL3.frequency(end));
    
            
	%% Get bottom depth
    
    %
    list_fields_aux = fieldnames(data_aux);
    %
    if any(strcmp(list_fields_aux, 'bottomdepth'))
        %
        str_h = 'bottomdepth';
	%
    elseif any(strcmp(list_fields_aux, 'bottomdepthfrompres'))
        %
        str_h = 'bottomdepthfrompres';
    else
        error('!!!')
    end
    
    %
    dataL3.bottomdepth(:, i) = data_aux.(str_h);
            
    
    %% Get Spp and zhab if not Spotter
    
    %
    if isfield(data_aux, 'Spp')
    
        %
        if strcmp(data_aux.instrument, "SoloD")
            %
            dataL3.zhab(i) = data_aux.zhab;
            
        %
        elseif strcmp(data_aux.instrument, "Signature") || ...
               strcmp(data_aux.instrument, "Aquadopp")

            %
            dataL3.zhab(i) = data_aux.zhab_transducer;
            dataL3.zhab_wavedir(i) = data_aux.zhab_vel;
        %
        else
            %
            error(['Unexpected instrument type when getting ' ...
                   'transducer height above the bottom!'])
        end

        %
        dataL3.Spp(:, l_getfreq, i) = data_aux.Spp(:, l_getfreq);
        
    end

    
	%% Get spectra of differenct quantities
    
    %
    dataL3.See(:, l_getfreq, i) = data_aux.See(:, l_getfreq);
    
    %
    dataL3.k(:, l_getfreq, i) = data_aux.k(:, l_getfreq);
    dataL3.cp(:, l_getfreq, i) = data_aux.cp(:, l_getfreq);
    dataL3.cg(:, l_getfreq, i) = data_aux.cg(:, l_getfreq);
    dataL3.E(:, l_getfreq, i) = data_aux.E(:, l_getfreq);
    dataL3.Ecg(:, l_getfreq, i) = data_aux.Ecg(:, l_getfreq);
    
    
    %% Get spectral quantities related to mean wave direction
    
    %
    if isfield(data_aux, list_vars_dir(1))
        %
        for i2 = 1:length(list_vars_dir)
            %
            dataL3.(list_vars_dir(i2))(:, l_getfreq, i) = ...
                                 data_aux.(list_vars_dir(i2))(:, l_getfreq);

        end
    end
    
    
 
end

   
%% Save structure

%
disp('Saving L3 merged data')

%
save(fullfile(dir_L3, 'roxsi_L3_allsites.mat'), "dataL3")
    
% print elapsed time
toc


