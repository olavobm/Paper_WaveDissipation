function [dataL4, pwfe] = pairwise_fe_getdatapair(dataL4, list_pairs)
%% [dataL4, pwfe] = PAIRWISE_FE_GETDATAPAIR(dataL4, list_pairs)
%
%   inputs
%       - dataL4: structure with L4 data.
%       - list_pairs: string array with location IDs for computing fe.
%
%   outputs
%       - dataL4: structure with L4 data (same as input).
%       - pwfe: structure variable with timeseries from instruments
%               pairs for computing fe.
%
%
% PAIRWISE_FE_GETDATAPAIR.m is a high-level function that gets data from
% a pair of instruments. These data is later used to compute fe between the
% instruments (in pairwise_fe_calculate.m).
%
%
% See also:
%   run_proc_pairwise_fe_allsites.m
%   pairwise_fe_calculate.m
%   find_ind_matchstring.m


%% Initialize structure and define parameters for RMS quantities

%
pwfe.Npairs = size(list_pairs, 1);
%
pwfe.locationID = list_pairs;

% Find indices matching locations between datasets
ind_match = find_ind_matchstring(dataL4.locationID, pwfe.locationID);

%
pwfe.ind_pairs = ind_match;


%% Add ChinaRock or Asilomar identifiers

%
pwfe.site = repmat("ChinaRock", size(list_pairs, 1), 1); 

%
for i = 1:size(list_pairs, 1)
	%
    id_aux = char(pwfe.locationID(i, 1));
    %
    if strcmp(id_aux(1), 'X')
        %
        pwfe.site(i, 1) = "Asilomar";
    end
end


%% Add variables to data structure

%
pwfe.freqlimsRMS = dataL4.seaswell.freqlims;

%
pwfe.dtime = dataL4.dtime;
pwfe.frequency = dataL4.frequency;


%%
% ----------------------------------------------
% ------ GET DATA FOR PAIRS OF LOCATIONS -------
% ----------------------------------------------
 
%
list_vars_spectra = ["See", "Ecg", "Fx", "Fy"];
list_vars_bulk_all = ["Hs", "Urms", "Ab"];

%
for i = 1:pwfe.Npairs
    
    %% Get location and water depth
    
    %
    pwfe.data(i).X = dataL4.X(:, pwfe.ind_pairs(i, :));
    pwfe.data(i).Y = dataL4.Y(:, pwfe.ind_pairs(i, :));

    
	%
    pwfe.data(i).bottomdepth = dataL4.bottomdepth(:, pwfe.ind_pairs(i, :));
    
    
    %% Get variables into pwfe structure
    
    %
    for i2 = 1:length(list_vars_spectra)
        
        pwfe.data(i).(list_vars_spectra(i2)) = ...
                dataL4.(list_vars_spectra(i2))(:, :, pwfe.ind_pairs(i, :));
        
    end
    
    %
    for i2 = 1:length(list_vars_bulk_all)
        %
        pwfe.data(i).(list_vars_bulk_all(i2)) = ...
                dataL4.seaswell.(list_vars_bulk_all(i2))(:, pwfe.ind_pairs(i, :));
    end
        
                         
	%% Average variables between sites
    
    %
    list_fields_aux = ["bottomdepth", "Hs", "Urms", "Ab"];
    
    %
    for i2 = 1:length(list_fields_aux)
        %
        pwfe.data(i).([char(list_fields_aux(i2)) '_avg']) = ...
                               mean(pwfe.data(i).(list_fields_aux(i2)), 2);
    end
    
    %% Get integrated flux
    
    %
    pwfe.data(i).fromEcg.Flux = dataL4.seaswell.Ecg(:, pwfe.ind_pairs(i, :));
    pwfe.data(i).fromFx.Flux = dataL4.seaswell.Fx(:, pwfe.ind_pairs(i, :));
    
    
    %% Compute cross-shore flux divergence
        
    %
    list_fields_aux = ["fromEcg", "fromFx"];

    %
    for i2 = 1:length(list_fields_aux)
        
        %
        diff_Flux = pwfe.data(i).(list_fields_aux(i2)).Flux(:, 2) - ...
                    pwfe.data(i).(list_fields_aux(i2)).Flux(:, 1);
        %
        diff_X = pwfe.data(i).X(:, 2) - pwfe.data(i).X(:, 1);
        
        %
        pwfe.data(i).(list_fields_aux(i2)).dFluxdx = diff_Flux ./ diff_X;

    end
                          

end
    
