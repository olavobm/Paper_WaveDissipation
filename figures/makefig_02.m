%% Make Figure 02 -- plot maps of sigma_h

clear
close all


%% Load bathymetry with roughness

%
dir_bathy = fullfile(paper_directory(), 'data', 'bathymetry');

% Using sigmah from pointcloud
bathymetry.Asilomar = load(fullfile(dir_bathy, 'bathy_longscale_Asilomar.mat'));
bathymetry.ChinaRock = load(fullfile(dir_bathy, 'bathy_longscale_ChinaRock.mat'));

%
bathymetry.Asilomar = bathymetry.Asilomar.bathymetry;
bathymetry.ChinaRock = bathymetry.ChinaRock.bathymetry;


%% Interpolate over gaps

%
for i_site = ["Asilomar", "ChinaRock"]
    
    %
    lok = ~isnan(bathymetry.(i_site).zmsl_stddev_deplaned);
    lgaps = ~isnan(lok);
    
    %
    [x_aux, y_aux] = meshgrid(bathymetry.(i_site).x, bathymetry.(i_site).y);
    
    %
    interp_aux = scatteredInterpolant(x_aux(lok), y_aux(lok), ...
                                      bathymetry.(i_site).zmsl_stddev_deplaned(lok));
    
	%
    z_aux = interp_aux(x_aux(lgaps), y_aux(lgaps));
    
    %
    bathymetry.(i_site).zmsl_stddev_deplaned(lgaps) = z_aux;
    
end


%% Load the 2m gridded elevation for land

%
bathymetry2m.ChinaRock = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_2m_ChinaRock.mat'));
bathymetry2m.ChinaRock = bathymetry2m.ChinaRock.bathymetry;
%
bathymetry2m.Asilomar = load(fullfile(paper_directory(), 'data', 'bathymetry', 'bathy_2m_Asilomar.mat'));
bathymetry2m.Asilomar = bathymetry2m.Asilomar.bathymetry;


%% Load instrument locations

%
mooringtable = load(fullfile(paper_directory(), 'proc', ...
                                                'deployment', ...
                                                'mooringtable_ROXSI2022.mat'));
mooringtable = mooringtable.mooringtable;

%
instrtable.Asilomar = mooringtable(strcmp(mooringtable.roxsiarray, "Asilomar"), :);
instrtable.ChinaRock = mooringtable(strcmp(mooringtable.roxsiarray, "ChinaRock"), :);


%% Define locations that will be plotted


%
pltsites.SoloD = ["A02"; "A04"; "A05"; "A06"; "A07"; "A08"; "A09"; ...
                  "B04"; "B06"; "B07"; "B09"; "B12"; "B14"; "B16"; "B18"; ...
                  "C02"; "C05"; "C06"; "C07"; "C08"; "C09"; ...
                  "D01"; "D02"; ...
                  "X02"; "X07"; "X08"; "X09"; "X10"; "X12"; "X14"];
%
pltsites.Spotter = ["B01"; "B03"; "B05"; "X01"; "X03"; "X04"];

%
pltsites.SmartMooring = ["E01"; "E02"; "E05"; "E07"; "E08"; "E09"; "E11"; "E13"];

%
pltsites.ADCP = ["B08"; "B10"; "B11"; "B13"; "B15"; "B17"; ...
                 "A01"; "C01"; ...
                 "E03"; "E04"; "E06"; "E12"; ...
                 "X05"; "X06"; "X11"; "X13"];
              
%
pltsites.allsites = [pltsites.SoloD; pltsites.Spotter; ...
                     pltsites.SmartMooring; pltsites.ADCP];
             

%% 
% ------------------------------------------------
% ------------ CREATE USEFUL COLORMAP ------------
% ------------------------------------------------

%%

%
figure, cmocean topo;
cmap_topo = get(gca, 'Colormap');
colorbar
close(gcf)


%%

%
zbottom = -22;
ztop = 10;

%
inds_trim_ocean = 24:124;
inds_trim_land = 129:256;


%%

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
cmap_new = [cmap_ocean; cmap_land];


%%
% -------------------------------------------------------------------
% --------------------------- MAKE FIGURE ---------------------------
% -------------------------------------------------------------------

%% Define color of symbols

%
clrcode.SoloD = [0, 0, 1];
clrcode.Spotter = [1, 1, 0];
clrcode.SmartMooring = [0, 1, 0];
clrcode.ADCP = [1, 0, 0];


%% Make figure


%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.4    0.2058    0.18    0.475];
%
haxs_ASL = axes('Position', [0.171, 0.565, 0.618, 0.618]);
haxs_CHR = axes('Position', [0.155, 0.125, 0.65, 0.65]);
%
hold(haxs_ASL, 'on')
hold(haxs_CHR, 'on')

%
mkSZout = 30;
mkSZin = 22;

% --------------------------------------------
% Plot pcolor of sigma_h

    %
    pcolor(haxs_ASL, bathymetry.Asilomar.x, bathymetry.Asilomar.y, bathymetry.Asilomar.zmsl_stddev_deplaned)
    shading(haxs_ASL, 'flat')
    
    %
    pcolor(haxs_CHR, bathymetry.ChinaRock.x, bathymetry.ChinaRock.y, bathymetry.ChinaRock.zmsl_stddev_deplaned)
    shading(haxs_CHR, 'flat')
    

% --------------------------------------------
% Plot instrument locations at Asilomar

    %
    for i = 1:size(instrtable.Asilomar, 1)
        %
        if any(strcmp(pltsites.allsites, instrtable.Asilomar.locationID(i)))
            %
            plot(haxs_ASL, instrtable.Asilomar.X(i), ...
                           instrtable.Asilomar.Y(i), '.k', 'MarkerSize', mkSZout)
            %
            if any(strcmp(pltsites.SoloD, instrtable.Asilomar.locationID(i)))
                clr_aux = clrcode.SoloD;
                
            %
            elseif any(strcmp(pltsites.Spotter, instrtable.Asilomar.locationID(i)))
                clr_aux = clrcode.Spotter;
                
            %
            elseif any(strcmp(pltsites.SmartMooring, instrtable.Asilomar.locationID(i)))
                clr_aux = clrcode.SmartMooring;
                
            %
            elseif any(strcmp(pltsites.ADCP, instrtable.Asilomar.locationID(i)))
                clr_aux = clrcode.ADCP;
                
            else
                error('!!!!')
            end
                
            %
            plot(haxs_ASL, instrtable.Asilomar.X(i), ...
                           instrtable.Asilomar.Y(i), ...
                           '.', 'Color', clr_aux, 'MarkerSize', mkSZin)
        end
        
    end
    


% --------------------------------------------
% Set colorbar and color limits

%
hcb = colorbar(haxs_CHR);
    %
    hcb.FontSize = 14;
    hcb.Position(1) = 0.81;
    hcb.Position(3) = 0.043;
    hcb.Position(4) = 0.83;
    %
    hcb.Label.Interpreter = 'Latex';
    hcb.Label.String = '$\sigma_h$ [m]';
    hcb.Label.FontSize = 20;
    
%
set([haxs_ASL, haxs_CHR], 'CLim', [0, 1.5])


% --------------------------------------------
% Set axes properties

%
set([haxs_ASL, haxs_CHR], 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')

%
set([haxs_ASL, haxs_CHR], 'DataAspectRatio', [1, 1, 1])
%
set(haxs_ASL, 'XLim', [-550, 0], 'YLim', [-100, 100])
set(haxs_CHR, 'XLim', [-550, 0], 'YLim', [-400, 400])
%
set(haxs_ASL, 'XTick', -400:200:0, 'YTick', [-100, 0, 100])
set(haxs_CHR, 'XTick', -400:200:0, 'YTick', -400:200:400)

%
set(haxs_ASL, 'YTick', [-50, 0, 50], 'XTickLabel', [])

%
hxlbl_CHR = xlabel(haxs_CHR, '$x$ [m]', 'Interpreter', 'Latex', 'FontSize', 18);
hylbl_CHR = ylabel(haxs_CHR, '$y$ [m]', 'Interpreter', 'Latex', 'FontSize', 18);
hylbl_ASL = ylabel(haxs_ASL, '$y$ [m]', 'Interpreter', 'Latex', 'FontSize', 18);
%
hxlbl_CHR.Position(2) = -465;
hylbl_CHR.Position(1) = -635;
hylbl_ASL.Position(1) = -635;

%
set([haxs_ASL, haxs_CHR], 'Color', 0.3.*[1, 1, 1])


% --------------------------------------------
% --------------------------------------------
% Make another axes to overlay the land

% ------------------------
%
haxs_CHR_land = axes('Position', [0.1, 0.1, 0.2, 0.2]);
hold(haxs_CHR_land, 'on')
haxs_CHR_land.Position = haxs_CHR.Position;
%
set(haxs_CHR_land, 'DataAspectRatio', [1, 1, 1])
set(haxs_CHR_land, 'XLim', [-550, 0], 'YLim', [-400, 400])
%
z_bla = bathymetry2m.ChinaRock.zmsl_mean;
z_bla(z_bla < 0) = NaN;
%
pcolor(haxs_CHR_land, bathymetry2m.ChinaRock.x, bathymetry2m.ChinaRock.y, z_bla)
shading flat
%
set(haxs_CHR_land, 'ColorMap', 0.35.*[1, 1, 1])


% ------------------------
%
haxs_ASL_land = axes('Position', [0.1, 0.1, 0.2, 0.2]);
hold(haxs_ASL_land, 'on')
haxs_ASL_land.Position = haxs_ASL.Position;
%
set(haxs_ASL_land, 'DataAspectRatio', [1, 1, 1])
set(haxs_ASL_land, 'XLim', [-550, 0], 'YLim', [-100, 100])
%
z_bla = bathymetry2m.Asilomar.zmsl_mean;
z_bla(z_bla < 0) = NaN;
%
pcolor(haxs_ASL_land, bathymetry2m.Asilomar.x, bathymetry2m.Asilomar.y, z_bla)
%
shading flat
%
set(haxs_ASL_land, 'ColorMap', 0.35.*[1, 1, 1])


% ------------------------
%
haxs_ASL_land.Visible = 'off';
haxs_CHR_land.Visible = 'off';


    % ------------------------
    % Plot instrument locations at China Rock

    %
    for i = 1:size(instrtable.ChinaRock, 1)
        %
        if any(strcmp(pltsites.allsites, instrtable.ChinaRock.locationID(i)))
            %
            plot(haxs_CHR_land, instrtable.ChinaRock.X(i), ...
                                instrtable.ChinaRock.Y(i), '.k', 'MarkerSize', mkSZout)
            %
            if any(strcmp(pltsites.SoloD, instrtable.ChinaRock.locationID(i)))
                clr_aux = clrcode.SoloD;

            %
            elseif any(strcmp(pltsites.Spotter, instrtable.ChinaRock.locationID(i)))
                clr_aux = clrcode.Spotter;

            %
            elseif any(strcmp(pltsites.SmartMooring, instrtable.ChinaRock.locationID(i)))
                clr_aux = clrcode.SmartMooring;

            %
            elseif any(strcmp(pltsites.ADCP, instrtable.ChinaRock.locationID(i)))
                clr_aux = clrcode.ADCP;

            else
                error('!!!!')
            end

            %
            plot(haxs_CHR_land, instrtable.ChinaRock.X(i), ...
                                instrtable.ChinaRock.Y(i), ...
                                '.', 'Color', clr_aux, 'MarkerSize', mkSZin)
        end
    end


% -----------------------------------------------
% Letter labeling

%
xrect = [-550, -550, -480, -480, -550];

%
fill(haxs_ASL, xrect, [30, 100, 100, 30, 30], 'w')
fill(haxs_CHR, xrect, [330, 400, 400, 330, 330], 'w')

%
text(haxs_ASL, -542, 65, 'a)', 'Interpreter', 'Latex', 'FontSize', 18)
text(haxs_CHR, -542, 365, 'b)', 'Interpreter', 'Latex', 'FontSize', 18)

%
hcb.Label.Position(1) = 2;


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure02.pdf'), 'Resolution', 300)



%%
return



%% Compute statistics of standard deviation

%
quantiles_range = [0.10, 0.90];


% ------------------------------
%
xlimsASL = [-550, 0];
ylimsASL = [-100, 100];

%
[xgridASL, ygridASL] = meshgrid(bathymetry.Asilomar.x, bathymetry.Asilomar.y);

%
lavgASL = (xgridASL >= xlimsASL(1)) & (xgridASL <= xlimsASL(2)) & ...
          (ygridASL >= ylimsASL(1)) & (ygridASL <= ylimsASL(2)) & ...
          (bathymetry.Asilomar.zmsl_mean <= 0) & ...
          ~isnan(bathymetry.Asilomar.zmsl_stddev_deplaned);
%
kbASL = bathymetry.Asilomar.zmsl_stddev_deplaned(lavgASL);

%
mean(kbASL)

%
quantile(kbASL, quantiles_range)

% ------------------------------
%
xlimsCHR = [-550, 0];
ylimsCHR = [-400, 400];

%
[xgridCHR, ygridCHR] = meshgrid(bathymetry.ChinaRock.x, bathymetry.ChinaRock.y);

%
lavgCHR = (xgridCHR >= xlimsCHR(1)) & (xgridCHR <= xlimsCHR(2)) & ...
          (ygridCHR >= ylimsCHR(1)) & (ygridCHR <= ylimsCHR(2)) & ...
          (bathymetry.ChinaRock.zmsl_mean <= 0) & ...
          ~isnan(bathymetry.ChinaRock.zmsl_stddev_deplaned);
%
kbCHR = bathymetry.ChinaRock.zmsl_stddev_deplaned(lavgCHR);

%
mean(kbCHR)

%
quantile(kbCHR, quantiles_range)




