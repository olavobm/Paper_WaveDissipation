%% Compute dh_bathy and adjust spectral estimates from pressure
% sensors with the effective depth correction (Marques et al., 2024).

close all
clear


%%
% -----------------------------------------------------
% ------ GET DEPTH CORRECTION AT INSTRUMENT SITES -----
% -----------------------------------------------------


%% Load L3 data

%
dir_dataL3 = fullfile(paper_directory(), 'data', 'level_3');

%
dataL3 = load(fullfile(dir_dataL3, 'roxsi_L3_allsites.mat'));
dataL3 = dataL3.dataL3;


%% Save file with a different name to keep results with
% observed depth for QC purposes.

%
file_hp_output = 'roxsi_L3_allsites_hp.mat';

% Check if file already exists
dir_info = dir(fullfile(dir_dataL3, file_hp_output));

% If it doesn't exist, save copy
if isempty(dir_info)
    %
    save(fullfile(dir_dataL3, file_hp_output), 'dataL3')

else
    % If it does exit, then terminate this script
    
    %
    file_aux = mfilename;
    %
    warning(['heff correction has already been applied. ' ...
             'Terminating ' file_aux '.m.'])
    return
end


%% Load bathymetry with depth corrections

%
dir_bathy = fullfile(paper_directory(), 'data', 'bathymetry');

%
bathymetry.Asilomar = load(fullfile(dir_bathy, 'bathy_dhbathy_Asilomar.mat'));
bathymetry.ChinaRock = load(fullfile(dir_bathy, 'bathy_dhbathy_ChinaRock.mat'));
%
bathymetry.Asilomar = bathymetry.Asilomar.bathymetry;
bathymetry.ChinaRock = bathymetry.ChinaRock.bathymetry;


%% Compute heff correction

%
dataL3.dh_bathy = NaN(length(dataL3.locationID), 1);

%
for i = 1:length(dataL3.locationID)
    
    %
	dh_aux = interp2(bathymetry.(dataL3.site(i)).x, ...
                     bathymetry.(dataL3.site(i)).y, ...
                     bathymetry.(dataL3.site(i)).dh_bathy, ...
                     mean(dataL3.X(:, i), 'omitnan'), ...
                     mean(dataL3.Y(:, i), 'omitnan'));
    
	%
    if isnan(dh_aux)
        dh_aux = 0;   % ~8 with NaN, but only really important for E08 and E09
    end
    
	%
    dataL3.dh_bathy(i) = dh_aux;
                            
end


%%
% ------------------------------------------------------------
% ------------------- CONVERT HEFF TO HOBS -------------------
% ------------------------------------------------------------

%
Nlocs = length(dataL3.locationID);


%% Get h and k with local depth

%
k_hp = dataL3.k;

%
hp_array = reshape(dataL3.bottomdepth, [length(dataL3.dtime), 1, length(dataL3.locationID)]);
hp_array = repmat(hp_array, [1, length(dataL3.frequency), 1]);


%% Recalculate bottom depth            

%
for i = 1:Nlocs
    %
	dataL3.bottomdepth(:, i) = dataL3.bottomdepth(:, i) + dataL3.dh_bathy(i);
end

%
heff_array = reshape(dataL3.bottomdepth, [length(dataL3.dtime), 1, length(dataL3.locationID)]);
heff_array = repmat(heff_array, [1, length(dataL3.frequency), 1]);


%% Recompute k

%
tic
%
for i1 = 1:Nlocs
    
    %
    for i2 = 1:length(dataL3.dtime)
        %
        dataL3.k(i2, :, i1) = wave_freqtok(dataL3.frequency, ...
                                           dataL3.bottomdepth(i2, i1));
        
    end
    
end
toc


%% Compute transfer function from See with hp and See with heff

%
zhab_hp = dataL3.zhab;
zhab_heff = zeros(size(zhab_hp));

%
zhab_hp_array = reshape(zhab_hp, [1, 1, Nlocs]);
zhab_heff_array = zeros([1, 1, Nlocs]);

%
lchange = (dataL3.dh_bathy > (-dataL3.zhab));
% % Old code, should only make a tiny difference
% lchange = (dataL3.dhbathy > 0);

%
for i = 1:Nlocs
    if lchange(i)
        %
        zhab_heff(i) = dataL3.dh_bathy(i) + dataL3.zhab(i);
        %
        zhab_heff_array(:, :, i) = zhab_heff(i);
    end
end

%
kh_term = cosh(dataL3.k .* heff_array) ./ cosh(k_hp .* hp_array);
kzhab_term = cosh(k_hp .* zhab_hp_array) ./ cosh(dataL3.k .* zhab_heff_array);

%
TF_Seeobs_to_Seeheff = (kh_term .* kzhab_term).^2;


% Don't apply to Spotters (and SmartMoorings)
ldontapply = strcmp(dataL3.instrument, "Spotter") | ...
             strcmp(dataL3.instrument, "SmartMooring");
% ldontapply should be the same vector as isnan(zhab_hp)
%
TF_Seeobs_to_Seeheff(:, :, ldontapply) = 1;

%
zhab_heff(ldontapply) = 0;

%
dataL3.zhab = zhab_heff;


%% Compute See for effective depth

%
dataL3.See = dataL3.See .* TF_Seeobs_to_Seeheff;


%% Recompute cp and cg

%
for i1 = 1:length(dataL3.locationID)    
    %
    for i2 = 1:length(dataL3.dtime)
        
        %
        dataL3.cp(i2, :, i1) = wave_cp(dataL3.k(i2, :, i1), ...
                                       dataL3.bottomdepth(i2, i1));
        
        %
        dataL3.cg(i2, :, i1) = wave_cg(dataL3.k(i2, :, i1), ...
                                       dataL3.bottomdepth(i2, i1));
        
    end
    
end


%% Recompute E and Ecg

%
dataL3.E = (1025 * 9.8) * dataL3.See;

%
dataL3.Ecg = (1025 * 9.8) * dataL3.cg .* dataL3.See;


%% Recompute Fx and Fy

%
dataL3.Fx = dataL3.Ecg .* dataL3.a1;
dataL3.Fy = dataL3.Ecg .* dataL3.b1;


%% Save (overwrite) L3 data file.

%
save(fullfile(dir_dataL3, 'roxsi_L3_allsites.mat'), 'dataL3')










