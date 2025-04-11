%% Make Figure 6 -- plot bathymetry and sigma_h between one instrument pair

clear
close all


%% Load bathymetry

%
bathymetry.ChinaRock = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_2m_ChinaRock.mat'));
bathymetry.ChinaRock = bathymetry.ChinaRock.bathymetry;
%
bathymetry.Asilomar = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_2m_Asilomar.mat'));
bathymetry.Asilomar = bathymetry.Asilomar.bathymetry;

%
bathymetrylongscale.ChinaRock = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_longscale_ChinaRock.mat'));
bathymetrylongscale.ChinaRock = bathymetrylongscale.ChinaRock.bathymetry;
%
bathymetrylongscale.Asilomar = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_longscale_Asilomar.mat'));
bathymetrylongscale.Asilomar = bathymetrylongscale.Asilomar.bathymetry;

% Using sigmah from pointcloud
bathymetry.ChinaRock.zmsl_stddev_deplaned = bathymetrylongscale.ChinaRock.zmsl_stddev_deplaned;
bathymetry.Asilomar.zmsl_stddev_deplaned = bathymetrylongscale.Asilomar.zmsl_stddev_deplaned;


%% Interpolate bathymetry over NaNs

%
list_aux = ["ChinaRock", "Asilomar"];

%
for i = 1:length(list_aux)
    
    %
    [xgrid, ygrid] = meshgrid(bathymetry.(list_aux(i)).x, bathymetry.(list_aux(i)).y);
    
    %
    lbad = isnan(bathymetry.(list_aux(i)).zmsl_mean);
    lgood = ~lbad;
    
    %
    interp_aux = scatteredInterpolant(xgrid(lgood), ygrid(lgood), ...
                                      bathymetry.(list_aux(i)).zmsl_mean(lgood));
                                  
	%
    bathymetry.(list_aux(i)).zmsl_mean(lbad) = interp_aux(xgrid(lbad), ygrid(lbad));
    
end


%% Load instrument locations

%
mooringtable = load(fullfile(paper_directory(), 'proc', ...
                                                'deployment', ...
                                                'mooringtable_ROXSI2022.mat'));
mooringtable = mooringtable.mooringtable;

%
instrtable.Asilomar = mooringtable(strcmp(mooringtable.roxsiarray, "Asilomar"), :);
instrtable.ChinaRock = mooringtable(strcmp(mooringtable.roxsiarray, "ChinaRock"), :);


%% Pick instrument pair

%
pick_pair = ["B11", "B12"];

%
ind_match_ID = find_ind_matchstring(instrtable.ChinaRock.locationID, pick_pair);


%% Get mean water at B11 and B12

%
dataL4 = load(fullfile(paper_directory(), 'data', 'level_4', 'roxsi_data_L4.mat'));
dataL4 = dataL4.dataL4;

%
ind_get_depth = find_ind_matchstring(dataL4.locationID, pick_pair);

%
depth_B11B12 = mean(dataL4.bottomdepth(:, ind_get_depth), 'omitnan');
%
X_B11B12 = mean(dataL4.X(:, ind_get_depth), 'omitnan');
Y_B11B12 = mean(dataL4.Y(:, ind_get_depth), 'omitnan');

%
depth_B11_mean = mean(depth_B11B12(:, 1), 'omitnan');
depth_B12_mean = mean(depth_B11B12(:, 2), 'omitnan');


%% Load data with sigmah per instrument pairs

%
dir_data = fullfile(paper_directory(), 'data', 'level_5');

%
load(fullfile(dir_data, 'dataL5_withQC.mat'))

%
ind_match_L5 = find_ind_matchstring(pwdissipation.locationID(:, 1), pick_pair(1));



%%
% -------------------------------------------
% ------------- CREATE COLORMAP -------------
% -------------------------------------------


%%

%
figure, cmocean topo;
cmap_topo = get(gca, 'Colormap');
colorbar
close(gcf)

%
zbottom = -10;
% ztop = -4;
ztop = -3;

%
% inds_trim_ocean = 24:128;
inds_trim_ocean = 24:124;
inds_trim_land = 129:256;

%
dz = 1;
%
Nlvls = length(zbottom:dz:ztop);
Nlvls_ocean = length(zbottom:dz:0);
Nlvls_land = length(0:dz:ztop);

%
inds_ocean_interp = linspace(inds_trim_ocean(1), inds_trim_ocean(end), Nlvls_ocean+1);
inds_land_interp = linspace(inds_trim_land(1), inds_trim_land(end), Nlvls_land);

%
cmap_ocean = interp1(inds_trim_ocean, cmap_topo(inds_trim_ocean, :), inds_ocean_interp);
cmap_land = interp1(inds_trim_land, cmap_topo(inds_trim_land, :), inds_land_interp);

%
cmap_new = cmap_ocean;



%%
% ------------------------------------------------------------
% ------------------------------------------------------------
% ----------------------- MAKE FIGURE ------------------------
% ------------------------------------------------------------
% ------------------------------------------------------------

%% Define limits for sigma_h

%
ylims_bathy = pwdissipation.data(ind_match_L5).yedges_bathy;

      
%% Get default Matlab color for instrument pair

Cmlab = get(groot,'DefaultAxesColorOrder');


%% Make figure

%
xWidth = 16;    % full page
yHeight = 13;

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.3777    0.1883    0.1845    0.3833];
%
hfig.PaperUnits = 'centimeters';
hfig.PaperPosition = [0, 0, xWidth, yHeight];

%
haxs_a = axes('Position', [0.175, 0.56, 0.8, 0.35]);
haxs_b = axes('Position', [0.175, 0.15, 0.8, 0.35]);
%
haxs_all = [haxs_a, haxs_b];
%
hold(haxs_all, 'on')

%
mkSZout = 38;
mkSZin = 28;

% ----------------------------------------------------
% Plot first panel with bathymetry

%
linylims = (bathymetry.ChinaRock.y >= ylims_bathy(1)) & ...
           (bathymetry.ChinaRock.y <= ylims_bathy(2));
ind_inybox = find(linylims);
%
ind_pick = round(linspace(2, length(ind_inybox)-1, 3));
ind_inyplt = ind_inybox(ind_pick);

            
    %
    plot(haxs_a, bathymetry.ChinaRock.x, ...
                 -mean(bathymetry.ChinaRock.zmsl_mean(linylims, :), 1), '-k', 'LineWidth', 2)
    x_aux = bathymetry.ChinaRock.x;
    y_aux = -mean(bathymetry.ChinaRock.zmsl_mean(linylims, :), 1);
    y_aux = y_aux(:);
    %
    x_vec_aux = [x_aux; x_aux(end); x_aux(1); x_aux(1)];
    y_vec_aux = [y_aux; 20; 20; y_aux(1)];
    
    %
	fill(haxs_a, x_vec_aux, y_vec_aux, 0.5.*[1, 1, 1])
    
    %
    plot(haxs_a, bathymetry.ChinaRock.x, ...
                 -mean(bathymetry.ChinaRock.zmsl_mean(linylims, :), 1) + ...
                   std(bathymetry.ChinaRock.zmsl_mean(linylims, :), 1), '--k', 'LineWidth', 2)
    %
    plot(haxs_a, bathymetry.ChinaRock.x, ...
                 -mean(bathymetry.ChinaRock.zmsl_mean(linylims, :), 1) + ...
                   -std(bathymetry.ChinaRock.zmsl_mean(linylims, :), 1), '--k', 'LineWidth', 2)


% ----------------------------------------------------
% Plot symbols and label B11, B12

%
plot(haxs_a, X_B11B12(1), depth_B11_mean, '.k', 'MarkerSize', 54)
plot(haxs_a, X_B11B12(2), depth_B12_mean, '.k', 'MarkerSize', 54)
%
plot(haxs_a, X_B11B12(1), depth_B11_mean, '.', 'Color', Cmlab(1, :), 'MarkerSize', 46)
plot(haxs_a, X_B11B12(2), depth_B12_mean, '.', 'Color', Cmlab(2, :), 'MarkerSize', 46)

%
text(haxs_a, X_B11B12(1)-14, 9.8, 'B11', 'Interpreter', 'Latex', 'FontSize', 20)
text(haxs_a, X_B11B12(2)-5, 8.5, 'B12', 'Interpreter', 'Latex', 'FontSize', 20)
    
                    
% ----------------------------------------------------
% Plot cross-shore lines of sigma_h

%
hp_lines = plot(haxs_b, bathymetry.ChinaRock.x, ...
                        bathymetry.ChinaRock.zmsl_stddev_deplaned(linylims, :).', ...
                        '-', 'LineWidth', 2);
	
%
for i = 1:length(hp_lines)
    %
    hp_lines(i).Color = [0.5.*[1, 1, 1], 0.7];
    
end
             
             
% ----------------------------------------------------
% Plot line showing where sigma_h is averaged
             
	%
    xlims_avg = pwdissipation.data(ind_match_L5).xedges_bathy;
    
             
	%
    linxavg = (bathymetry.ChinaRock.x >= xlims_avg(1)) & ...
              (bathymetry.ChinaRock.x <= xlims_avg(2));
	%
    sigmah_sub = bathymetry.ChinaRock.zmsl_stddev_deplaned(linylims, linxavg);
    sigmah_sub = sigmah_sub(:);
    %
    sigmah_avg = mean(sigmah_sub);
    
	%
    plot(haxs_b, xlims_avg, sigmah_avg.*[1, 1], '-k', 'LineWidth', 4)
            
    
    % % % Check these match
    % % [pwdissipation.data(ind_match_L5).zmsl_stddev_avg, sigmah_avg];
    
    %
    text(haxs_b, -244, 0.8, ['$\langle \sigma_h \rangle = ' num2str(sigmah_avg, '%.1f') '~\mathrm{m}$'], ...
                               'Interpreter', 'Latex', 'FontSize', 20)
    
% ----------------------------------------------------
% Set properties of axes

%
set(haxs_all, 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')

%
set(haxs_a, 'YLim', [3.5, 10.4])
set(haxs_b, 'YLim', [0.5, 1.6])

%
set(haxs_all, 'XLim', [-270, -195])
%
set(haxs_a, 'YDir', 'reverse')

%
set(haxs_all(1:end-1), 'XTickLabel', [])


% ----------------------------------------------------
% Add x/y labels

%
hlby_A = ylabel(haxs_a, '$\overline{h}$ [m]', 'Interpreter', 'Latex', 'FontSize', 16);
hlby_B = ylabel(haxs_b, '$\sigma_h$ [m]', 'Interpreter', 'Latex', 'FontSize', 16);
%
hlby_A.Position(1) = -276.5;
hlby_B.Position(1) = -276.5;
%
xlabel(haxs_b, '$x$ [m]', 'Interpreter', 'Latex', 'FontSize', 16)


% ----------------------------------------------------
% Letter labelling

%
xposlettter = -269;
%
xrect = [-270, -270, -262.5, -262.5, -270];
%
hf_a = fill(haxs_a, xrect, [3.5, 5.1, 5.1, 3.5, 3.5], 'w');
hf_b = fill(haxs_b, xrect, [1.28, 1.5, 1.5, 1.28, 1.28]+0.1, 'w');
%
text(haxs_a, xposlettter, 4.4, 'a)', 'Interpreter', 'Latex', 'FontSize', 20)
text(haxs_b, xposlettter, 1.48, 'b)', 'Interpreter', 'Latex', 'FontSize', 20)


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure06.pdf'), 'Resolution', 300)

                          
                          
                          
                          
                          
                          
                          
                          