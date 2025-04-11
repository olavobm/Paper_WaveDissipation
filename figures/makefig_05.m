%% Make Figure 5 -- timeseries showing fe calculation.

clear
close all


%% Load L5 data with quality control

%
load(fullfile(paper_directory(), 'data', 'level_5', 'dataL5_withQC.mat'))


%% Pick a pair of instrument sites

%
pick_pair = ["B11", "B12"];

%
ind_match_A = find_ind_matchstring(pwdissipation.locationID(:, 1), pick_pair(1));
ind_match_B = find_ind_matchstring(pwdissipation.locationID(:, 2), pick_pair(2));

%
if ind_match_A==ind_match_B
    %
    ind_match = ind_match_A;
else
    error('!!!')
end

%
str_1 = char(pick_pair(1));
str_2 = char(pick_pair(2));


%% Load bathymetry

%
bathymetry.ChinaRock = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_2m_ChinaRock.mat'));
bathymetry.ChinaRock = bathymetry.ChinaRock.bathymetry;
%
bathymetry.Asilomar = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_2m_Asilomar.mat'));
bathymetry.Asilomar = bathymetry.Asilomar.bathymetry;


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


%% Create colormap

%
figure, cmocean topo;
cmap_topo = get(gca, 'Colormap');
colorbar
close(gcf)

%
zbottom = -10;
% ztop = -4;
ztop = -1;

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
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------- MAKE FIGURE ----------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------


%%

%
pick_pair = ["B11", "B12"];


%
plt_sites = ["B11", "B12", "B12", "B13", "B14", "B15"];

%
ind_match_ID = find_ind_matchstring(instrtable.ChinaRock.locationID, plt_sites);
%
ind_match_ID_pair = find_ind_matchstring(instrtable.ChinaRock.locationID, pick_pair);


%% Get standard Matlab colors for coloring the instrument pair

Cmlab = get(groot, 'DefaultAxesColorOrder');


%% Define limits of x y axes

%
xlims_avg = [(instrtable.ChinaRock.X(ind_match_ID_pair(1))-10), ...
             (instrtable.ChinaRock.X(ind_match_ID_pair(2))+10)];
%
ylims_bathy = [-22, 4];


%% Make figure


%
dtime_lims = [datetime(2022, 06, 21, 18, 00, 00), ...
              datetime(2022, 07, 21, 04, 00, 00)];
dtime_lims.TimeZone = 'America/Los_Angeles';

%
h_avg = mean(pwdissipation.data(ind_match).bottomdepth, 1, 'omitnan');

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.3545    0.0525    0.25    0.7];
%
haxs_map = makeSubPlots(0.19, 0.15, 0.1, ...
                        0.03, 0.77, 0.02, 1, 1);
%
haxs = makeSubPlots(0.19, 0.15, 0.1, ...
                    0.25, 0.08, 0.0125, 1, 4);
%
hold(haxs_map, 'on')
hold(haxs, 'on')

    % ----------------------------------------------------
    % pcolor bathymetry
    
    %
    pcolor(haxs_map, bathymetry.ChinaRock.x, bathymetry.ChinaRock.y, bathymetry.ChinaRock.zmsl_mean)
    
%
shading(haxs_map, 'flat')

    % ----------------------------------------------------       
    % Plot rectangle for sigma_h
    
%
plot(haxs_map, pwdissipation.data(ind_match).xedges_bathy(1).*[1, 1], pwdissipation.data(ind_match).yedges_bathy, '-k', 'LineWidth', 1.5)
plot(haxs_map, pwdissipation.data(ind_match).xedges_bathy(2).*[1, 1], pwdissipation.data(ind_match).yedges_bathy, '-k', 'LineWidth', 1.5)
%
plot(haxs_map, pwdissipation.data(ind_match).xedges_bathy, pwdissipation.data(ind_match).yedges_bathy(1).*[1, 1], '-k', 'LineWidth', 1.5)
plot(haxs_map, pwdissipation.data(ind_match).xedges_bathy, pwdissipation.data(ind_match).yedges_bathy(2).*[1, 1], '-k', 'LineWidth', 1.5)


    % ----------------------------------------------------   
    % Plot instrument sites
    
    %
    plot(haxs_map, instrtable.ChinaRock.X(ind_match_ID(1:2)), ...
                   instrtable.ChinaRock.Y(ind_match_ID(1:2)), ...
                   '.', 'Color', 0.3.*[1, 1, 1], 'MarkerSize', 72)
    %
    plot(haxs_map, instrtable.ChinaRock.X(ind_match_ID(3:end)), ...
                   instrtable.ChinaRock.Y(ind_match_ID(3:end)), ...
                   '.', 'Color', 0.3.*[1, 1, 1], 'MarkerSize', 52)
               
    %
    plot(haxs_map, instrtable.ChinaRock.X(ind_match_ID_pair(1)), ...
                   instrtable.ChinaRock.Y(ind_match_ID_pair(1)), ...
                   '.', 'MarkerSize', 62, 'Color', Cmlab(1, :))
    %
    plot(haxs_map, instrtable.ChinaRock.X(ind_match_ID_pair(2)), ...
                   instrtable.ChinaRock.Y(ind_match_ID_pair(2)), ...
                   '.', 'MarkerSize', 62, 'Color', Cmlab(2, :))


    % -----------------------------------------------
    % Plot timeseries
    
    %
    plot(haxs(1), pwdissipation.dtime, squeeze(pwdissipation.data(ind_match).Hs), 'LineWidth', 2);
    plot(haxs(2), pwdissipation.dtime, squeeze(pwdissipation.data(ind_match).fromEcg.Flux)./1000, 'LineWidth', 2)
    %
    axes(haxs(3))
    yyaxis left
    plot(haxs(3), pwdissipation.dtime, -squeeze(pwdissipation.data(ind_match).fromEcg.dFluxdx), '-k', 'LineWidth', 2)
    yyaxis right
    hp_Urms = plot(haxs(3), pwdissipation.dtime, squeeze(pwdissipation.data(ind_match).Urms_avg).^3, '-r', 'LineWidth', 2);

    %
    plot(haxs(4), pwdissipation.dtime, pwdissipation.data(ind_match).fromEcg.fe_timeseries, '-k', 'LineWidth', 2);
    %
    plot(haxs(4), dtime_lims, pwdissipation.data(ind_match).fromEcg.fctefit.*[1, 1], '--r', 'LineWidth', 2.5);

    %
    text(haxs(4), pwdissipation.dtime(525), 1.5, '$\tilde{f}_e$', 'Color', 'r', 'Interpreter', 'Latex', 'FontSize', 19)

     
    % -----------------------------------------------
    % Make legend
    
    %
    hleg = legend(haxs(2), [char(pick_pair(1)) ', $\overline{h} = ' num2str(h_avg(1), '%.1f') '~\mathrm{m}$'], ...
                           [char(pick_pair(2)) ', $\overline{h} = ' num2str(h_avg(2), '%.1f') '~\mathrm{m}$'], ...
                           'Interpreter', 'Latex');
        hleg.FontSize = 16;
        hleg.Position = [0.3107    0.5008    0.3453    0.0697];

    
% -----------------------------------------------
% Set axes

%
set(haxs, 'FontSize', 13, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
set(haxs, 'XLim', pwdissipation.dtime([1, end]))
set(haxs(1:end-1), 'XTickLabel', [])
%
set(haxs, 'XLim', dtime_lims)
%
set(haxs(1), 'YLim', [0, 2])
set(haxs(2), 'YLim', [0, 17])
set(haxs(4), 'YLim', [0, 13])
%
set(haxs(1), 'YTick', [0, 0.5, 1, 1.5, 2])
%
both_yaxes = get(haxs(3), 'YAxis');
both_yaxes(1).Color = 'k';
both_yaxes(2).Color = 'r';
%
both_yaxes(1).Limits = [0, 125];
both_yaxes(2).Limits = [0, 0.05];

%
both_yaxes(1).TickValues = [0, 25, 50, 75, 100];
both_yaxes(2).TickValues = [0, 0.01, 0.02, 0.03, 0.04, 0.05];


% --------------------------
%
set(haxs_map, 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')

%
set(haxs_map, 'DataAspectRatio', [1, 1, 1])

%
set(haxs_map, 'XLim', [-280, -80], 'YLim', [-30, 30])


% ----------------------------------------------------
% Set colormap and colorbar

%
set(haxs_map, 'Colormap', cmap_new, 'CLim', [zbottom, ztop])

%
hcb = colorbar(haxs_map);
    hcb.Position(1) = 0.86;
    hcb.Position(2) = 0.8075;
    hcb.Position(3) = 0.025;
    hcb.Position(4) = 0.125;
    %
    hcb.Label.Interpreter = 'Latex';
    hcb.Label.String = '$\overline{h}$ [m]';
    hcb.Label.FontSize = 18;
    %
    hcb.Ticks = [-10, -8, -6, -4, -2];
    hcb.TickLabels = {'10', '8', '6', '4', '2'};
    

% -----------------------------------------------
% Add labels

%
lblFS = 18;
%
ylabel(haxs(1), '$H_s$ [m]', 'Interpreter', 'Latex', 'FontSize', lblFS)
ylabel(haxs(2), '$F$ [kW m$^{-1}$]', 'Interpreter', 'Latex', 'FontSize', lblFS)
ylabel(haxs(4), '$f_e$', 'Interpreter', 'Latex', 'FontSize', lblFS)

%
both_yaxes(1).Label.Interpreter = 'Latex';
both_yaxes(2).Label.Interpreter = 'Latex';
%
both_yaxes(1).Label.String = {'$-\mathrm{d} F/ \mathrm{d} x$';'[W m$^{-2}$]'};
both_yaxes(2).Label.String = '$\langle U_{\mathrm{rms}} \rangle^3$ [m$^3$ s$^{-3}$]';
%
both_yaxes(1).Label.FontSize = lblFS;
both_yaxes(2).Label.FontSize = lblFS;


%
hxlb = xlabel(haxs_map, '$x$ [m]', 'Interpreter', 'Latex', 'FontSize', 16);
ylabel(haxs_map, '$y$ [m]', 'Interpreter', 'Latex', 'FontSize', 16)
%
hxlb.Position(2) = hxlb.Position(2) + 8;


% -----------------------------------------------
% Annotation for instrument sites

text(haxs_map(1), -260, 20, 'B11', 'Interpreter', 'Latex', 'FontSize', 18)
text(haxs_map(1), -225, 20, 'B12', 'Interpreter', 'Latex', 'FontSize', 18)
text(haxs_map(1), -198, -3, 'B13', 'Interpreter', 'Latex', 'FontSize', 18)
text(haxs_map(1), -160, 5, 'B14', 'Interpreter', 'Latex', 'FontSize', 18)
text(haxs_map(1), -110, -15, 'B15', 'Interpreter', 'Latex', 'FontSize', 18)


% -----------------------------------------------
% Letter labelling

%
xtxt_1 = pwdissipation.dtime(155);
xtxt_2 = pwdissipation.dtime(184);

%
fill(haxs_map(1), [-280, -280, -262, -262, -280], [10, 30, 30, 10, 10], 'w')
%
text(haxs_map(1), -278, 20, 'a)', 'FontSize', 18)

%
text(haxs(1), xtxt_1, 1.825, 'b)', 'FontSize', 18)
text(haxs(2), xtxt_1, 15.5, 'c)', 'FontSize', 18)
text(haxs(3), xtxt_1, 0.045, 'd)', 'FontSize', 18)
text(haxs(4), xtxt_1, 11.7, 'e)', 'FontSize', 18)


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure05.pdf'), 'Resolution', 300)


