%% Make Figure 1 -- maps

clear
close all


%% Load bathymetry

%
bathymetry.ChinaRock = load(fullfile(paper_directory(), ...
                                     'data', 'bathymetry', ...
                                     'bathy_2m_ChinaRock.mat'));
bathymetry.ChinaRock = bathymetry.ChinaRock.bathymetry;
%
bathymetry.Asilomar = load(fullfile(paper_directory(), ...
                                    'data', 'bathymetry', ...
                                    'bathy_2m_Asilomar.mat'));
bathymetry.Asilomar = bathymetry.Asilomar.bathymetry;


%% Interpolate over NaNs

%
for i_site = ["Asilomar", "ChinaRock"]
    
    %
    lok = ~isnan(bathymetry.(i_site).zmsl_stddev_deplaned);
    lgaps = ~isnan(lok);
    
    %
    [x_aux, y_aux] = meshgrid(bathymetry.(i_site).x, bathymetry.(i_site).y);
    
    %
    interp_aux = scatteredInterpolant(x_aux(lok), y_aux(lok), ...
                                      bathymetry.(i_site).zmsl_mean(lok));
    
	%
    z_aux = interp_aux(x_aux(lgaps), y_aux(lgaps));
    
    %
    bathymetry.(i_site).zmsl_mean(lgaps) = z_aux;
    
end


%% Load photos

%
photo_MP_georect = imread(fullfile(paper_directory(), 'figures', 'MP_googleearth_highres.jpg'));
%
photo_ChinaRock = imread(fullfile(paper_directory(), 'figures', 'photo_shoreline.jpg'));


%% Load instrument locations

%
mooringtable = load(fullfile(paper_directory(), 'proc', ...
                                                'deployment', ...
                                                'mooringtable_ROXSI2022.mat'));
mooringtable = mooringtable.mooringtable;

%
instrtable.Asilomar = mooringtable(strcmp(mooringtable.roxsiarray, "Asilomar"), :);
instrtable.ChinaRock = mooringtable(strcmp(mooringtable.roxsiarray, "ChinaRock"), :);


%% Set output directory for figures

dir_output = fullfile(paper_directory(), 'figures');


%% Load georectification parameters for image with the Monterey Peninsula

%
MPgeorect = load(fullfile(paper_directory(), 'figures', 'utils', 'georect_MP.mat'));
MPgeorect = MPgeorect.coefsgeorect;


%% 
% ------------------------------------------------
% --- SET SITE LOCATIONS PER GROUP OF DATASET ----
% ------------------------------------------------

%%

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


% the transition in colors is in the middle --  indices 128-129.
% first the ocean colors, then the land


%% Trim/interpolate to create colormap just for the ocean

%
zbottom = -24;
ztop = 0;

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



%%
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -------------------------- MAKE FIGURE --------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------


%%

%
clrcode.SoloD = [0, 0, 1];
clrcode.Spotter = [1, 1, 0];
clrcode.SmartMooring = [0.85, 0.325, 0.098];
clrcode.SmartMooring = [0, 1, 0];
clrcode.ADCP = [1, 0, 0];


%% Latitude/longitude limits for large-scale map

% West coast
lat_lims = [22, 54];
lon_lims = [-130, -108];


%% x/y limits for smaller scale maps

%
xlims_ASL_CHR = [-860, 20];
ylims_ASL = [-100, +100];
ylims_CHR = [-450, +450];

%
plthlID = ["B03", "B13"];

%
% % xWidth = 8;
xWidth = 16;    % full page
yHeight = 13;

%
mkSZout = 38;
mkSZin = 28;


%% Make figure

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.316363636363636   0.180000000000000   0.352272727272727   0.516666666666667];
% % %
% % hfig.PaperUnits = 'centimeters';
% % hfig.PaperPosition = [0, 0, xWidth, yHeight];

%
haxs_a = axes('Position', [0, 0.385, 0.495 0.495]);   % Google Earth/MP
haxs_b = axes('Position', [0.26, 0.41, 0.175, 0.175]);   % m_map / USA
haxs_c = axes('Position', [0.0735, 0.0375, 0.348, 0.348]);   % ChinaRock photo
%
haxs_d = axes('Position', [0.515, 0.545, 0.468, 0.468]);   % Asilomar map
haxs_e = axes('Position', [0.442, 0.09, 0.6135, 0.6135]);   % ChinaRock map
%
hold(haxs_b, 'on')
hold(haxs_d, 'on')
hold(haxs_e, 'on')

    % ----------------------------------------------------
    % Google Earth with MP
    
    %
    image(haxs_a, photo_MP_georect)
    %
    hold(haxs_a, 'on')
    %
    set(haxs_a, 'DataAspectRatio', [1, 1, 1])
    
    %
    set(haxs_a, 'XGrid', 'on', 'YGrid', 'on')
    %
    set(haxs_a, 'XLim', [3300, 5000], 'YLim', [1600, 3500])

    
    % % % --------------------
    % % % Plot instrument locations on Google map panel
    % % for i = 1:length(pltsites.allsites)
    % %     %
    % %     if any(strcmp(instrtable.Asilomar.locationID, pltsites.allsites(i)))
    % %         %
    % %         siteplt = 'Asilomar';
    % %     else
    % %         siteplt = 'ChinaRock';
    % %     end
    % % 
    % %     %
    % %     imatch_aux = find(strcmp(instrtable.(siteplt).locationID, pltsites.allsites(i)));
    % % 
    % %     %
    % %     x_aux = MPgeorect.fcn_x(instrtable.(siteplt).longitude(imatch_aux));
    % %     y_aux = MPgeorect.fcn_y(instrtable.(siteplt).latitude(imatch_aux));
    % %     %
    % %     plot(haxs_a, x_aux, y_aux, '.k', 'MarkerSize', 4)
    % % 
    % % 
    % % end

    
    % --------------------
    % Plot rectangles in Google image figure
    
    %
    [latpt1, lonpt1] = ROXSI_xytoll('ChinaRock', xlims_ASL_CHR(1), ylims_CHR(1));
    [latpt2, lonpt2] = ROXSI_xytoll('ChinaRock', xlims_ASL_CHR(1), ylims_CHR(2));
    %
    [latpt3, lonpt3] = ROXSI_xytoll('ChinaRock', xlims_ASL_CHR(2), ylims_CHR(2));
    [latpt4, lonpt4] = ROXSI_xytoll('ChinaRock', xlims_ASL_CHR(2), ylims_CHR(1));
    %
    rect.CHR.latitude = [latpt1, latpt2, latpt3, latpt4, latpt1];
    rect.CHR.longitude = [lonpt1, lonpt2, lonpt3, lonpt4, lonpt1];
    
    %
    rect.CHR.x = MPgeorect.fcn_x(rect.CHR.longitude);
    rect.CHR.y = MPgeorect.fcn_y(rect.CHR.latitude);
    
    %
    [latpt1, lonpt1] = ROXSI_xytoll('Asilomar', xlims_ASL_CHR(1), ylims_ASL(1));
    [latpt2, lonpt2] = ROXSI_xytoll('Asilomar', xlims_ASL_CHR(1), ylims_ASL(2));
    %
    [latpt3, lonpt3] = ROXSI_xytoll('Asilomar', xlims_ASL_CHR(2), ylims_ASL(2));
    [latpt4, lonpt4] = ROXSI_xytoll('Asilomar', xlims_ASL_CHR(2), ylims_ASL(1));
    %
    rect.ASL.latitude = [latpt1, latpt2, latpt3, latpt4, latpt1];
    rect.ASL.longitude = [lonpt1, lonpt2, lonpt3, lonpt4, lonpt1];
    
    %
    rect.ASL.x = MPgeorect.fcn_x(rect.ASL.longitude);
    rect.ASL.y = MPgeorect.fcn_y(rect.ASL.latitude);

    %
    lnWdaux = 2;
    %
    plot(haxs_a, rect.CHR.x(1:2), rect.CHR.y(1:2), '-r', 'LineWidth', lnWdaux)
    plot(haxs_a, rect.CHR.x(2:3), rect.CHR.y(2:3), '-r', 'LineWidth', lnWdaux)
    plot(haxs_a, rect.CHR.x(3:4), rect.CHR.y(3:4), '-r', 'LineWidth', lnWdaux)
    plot(haxs_a, rect.CHR.x([1, 4]), rect.CHR.y([1, 4]), '-r', 'LineWidth', lnWdaux)
    %
    plot(haxs_a, rect.ASL.x(1:2), rect.ASL.y(1:2), '-r', 'LineWidth', lnWdaux)
    plot(haxs_a, rect.ASL.x(2:3), rect.ASL.y(2:3), '-r', 'LineWidth', lnWdaux)
    plot(haxs_a, rect.ASL.x(3:4), rect.ASL.y(3:4), '-r', 'LineWidth', lnWdaux)
    plot(haxs_a, rect.ASL.x([1, 4]), rect.ASL.y([1, 4]), '-r', 'LineWidth', lnWdaux)
    
    
    % --------------------
    % Set longitude tick marks
    %
    lon_min = -58 : 4: -54;
    lon_ticks = -121 + (lon_min./60);
    %
    xticks_MP = MPgeorect.fcn_x(lon_ticks);
    %
    xticklabel_MP = cell(size(xticks_MP));
    %
    for i = 1:length(xticklabel_MP)
        xticklabel_MP{i} = ['121$^\circ$ ' num2str(abs((lon_ticks(i) + 121).*60), '%.0f') '''W'];
    end
    
    %
    set(haxs_a, 'XTick', xticks_MP, 'XTickLabel', xticklabel_MP)
    
    
    % --------------------
    % Set latitude tick marks (with stacked symbols)
    %
    lat_min = 28 : 2 : 40;
    lat_ticks = 36 + (lat_min./60);
    %
    yticks_MP = MPgeorect.fcn_y(lat_ticks);
    yticks_MP = fliplr(yticks_MP);
    
    %
    yticklabel_MP_row_1 = cell(size(yticks_MP));
    yticklabel_MP_row_2 = cell(size(yticks_MP));
    
    %
    for i = length(yticklabel_MP_row_1):-1:1
        yticklabel_MP_row_1{i} = '36^\circ \hspace{-0.28cm}';
        yticklabel_MP_row_2{i} = ['$$' num2str((lat_ticks(i) - 36).*60, '%.0f') '''N \hspace{-0.4cm} $$ '];
    end
    yticklabel_MP_row_1 = fliplr(yticklabel_MP_row_1);
    yticklabel_MP_row_2 = fliplr(yticklabel_MP_row_2);
    
    %
    yticklabel_MP_row = cell(size(yticks_MP));
    %
    for i = 1:length(yticklabel_MP_row)
        yticklabel_MP_row{i} = ['$$\begin{array}{c} ' yticklabel_MP_row_1{i} ' \\ ' yticklabel_MP_row_2{i} '  \\ \end{array}$$'];
    end
    
    %
    set(haxs_a, 'YTick', yticks_MP, 'YTickLabel', yticklabel_MP_row)

    %
    set(haxs_a, 'TickLabelInterpreter', 'latex')
    
    
    % ----------------------------------------------------
    % m_map with West coast
    
    %
    axes(haxs_b)

        %
        m_proj('robinson', 'lat', lat_lims, 'lon', lon_lims)
        m_grid('box', 'on', 'xticklabels', [], 'yticklabels', [])

        %
        fill_clr = 0.35.*[1, 1, 1];
        %
        m_coast('patch', fill_clr);
        
        %
        lat_ref = 36 + 35.5/60;
        lon_ref = -121 - 56/60;
        
        %
        m_plot(lon_ref, lat_ref, 'pr', 'MarkerFaceColor', 'r', 'MarkerSize', 18)
        
	% ----------------------------------------------------
    % China Rock photo
        
    %
    image(haxs_c, photo_ChinaRock)
    %
    hold(haxs_c, 'on')
    
    %
    set(haxs_c, 'DataAspectRatio', [1, 1, 1])
    
    %
    set(haxs_c, 'XTickLabel', [], 'YTickLabel', [])

    %
    xlim(haxs_c, [0.5, 1500])
    ylim(haxs_c, [100, 913.5])
    
    % ----------------------------------------------------
    % Asilomar
    
    %
    z_bla = bathymetry.Asilomar.zmsl_mean;
    z_bla(z_bla >= 0) = NaN;
    
    %
    pcolor(haxs_d, bathymetry.Asilomar.x, bathymetry.Asilomar.y, z_bla)
    %
    shading(haxs_d, 'flat')
    

    % ----------------------------------------------------
    % China Rock
    
    %
    z_bla = bathymetry.ChinaRock.zmsl_mean;
    z_bla(z_bla >= 0) = NaN;
    
    %
    pcolor(haxs_e, bathymetry.ChinaRock.x, bathymetry.ChinaRock.y, z_bla)
    %
    shading(haxs_e, 'flat')
    
    
% ----------------------------------------------------
% Set colormaps and colorbar


% -----------------
%
set([haxs_d, haxs_e], 'Colormap', cmap_ocean, 'CLim', [zbottom, ztop])


% -----------------
%
xcb = 0.55;
wdb = 0.03;
cbtcksFS = 13;
    
%
hcb_1 = colorbar(haxs_e, 'NorthOutside');
    %
    hcb_1.FontSize = cbtcksFS;
	%
    hcb_tick_marks = hcb_1.Ticks;
    hcb_tick_marks = abs(hcb_tick_marks);
    %
    for i = 1:length(hcb_1.TickLabels)
        hcb_1.TickLabels{i} = num2str(hcb_tick_marks(i));
    end
    
    
    %
    hcb_1.Position(1) = haxs_d.Position(1);
    hcb_1.Position(3) = haxs_d.Position(3);
    hcb_1.Position(2) = 0.86;
    hcb_1.Position(4) = 0.02;
    
	%
    hcb_1.Label.Interpreter = 'Latex';
    hcb_1.Label.String = '$\overline{h}$ [m]';
    hcb_1.Label.FontSize = 16;
    hcb_1.Label.Rotation = 0;
    
    
% ----------------------------------------------------
% Set axes propertiers for small-scale maps

%
set([haxs_d, haxs_e], 'DataAspectRatio', [1, 1, 1])

%
set(haxs_d, 'XLim', xlims_ASL_CHR, 'YLim', ylims_ASL)
%
set(haxs_e, 'XLim', xlims_ASL_CHR, 'YLim', ylims_CHR)

%
set([haxs_d, haxs_e], 'FontSize', 12)

%
set(haxs_d, 'XTickLabel', [])


% ----------------------------------------------------
% Set other axes properties

%
set([haxs_b, haxs_c], 'FontSize', 14)
%
set([haxs_b, haxs_c], 'Color', 0.*[1, 1, 1])


% ----------------------------------------------------
%
hylbl_c = ylabel(haxs_d, '$y$ [m]', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel(haxs_e, '$x$ [m]', 'Interpreter', 'Latex', 'FontSize', 16)
hylbl_d = ylabel(haxs_e, '$y$ [m]', 'Interpreter', 'Latex', 'FontSize', 16);
%
hylbl_c.Position(1) = -930;
hylbl_d.Position(1) = -930;

    
% --------------------------------------------
% Make another axes to overlay the land

% -------------
%
haxs_ASL_land = axes('Position', [0.1, 0.1, 0.2, 0.2]);
hold(haxs_ASL_land, 'on')
haxs_ASL_land.Position = haxs_d.Position;
%
set(haxs_ASL_land, 'DataAspectRatio', [1, 1, 1])
set(haxs_ASL_land, 'XLim', xlims_ASL_CHR, 'YLim', ylims_ASL)
%
z_bla = bathymetry.Asilomar.zmsl_mean;
z_bla(z_bla < 0) = NaN;
%
pcolor(haxs_ASL_land, bathymetry.Asilomar.x, bathymetry.Asilomar.y, z_bla)
shading flat
%
set(haxs_ASL_land, 'ColorMap', 0.35.*[1, 1, 1])


% -------------
%
haxs_CHR_land = axes('Position', [0.1, 0.1, 0.2, 0.2]);
hold(haxs_CHR_land, 'on')
haxs_CHR_land.Position = haxs_e.Position;
%
set(haxs_CHR_land, 'DataAspectRatio', [1, 1, 1])
set(haxs_CHR_land, 'XLim', xlims_ASL_CHR, 'YLim', ylims_CHR)
%
z_bla = bathymetry.ChinaRock.zmsl_mean;
z_bla(z_bla < 0) = NaN;
%
pcolor(haxs_CHR_land, bathymetry.ChinaRock.x, bathymetry.ChinaRock.y, z_bla)
shading flat
%
set(haxs_CHR_land, 'ColorMap', 0.35.*[1, 1, 1])


% -------------
    %
    for i = 1:size(instrtable.Asilomar, 1)
        %
        if any(strcmp(pltsites.allsites, instrtable.Asilomar.locationID(i)))
            %
            plot(haxs_ASL_land, instrtable.Asilomar.X(i), ...
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
            plot(haxs_ASL_land, instrtable.Asilomar.X(i), ...
                         instrtable.Asilomar.Y(i), ...
                         '.', 'Color', clr_aux, 'MarkerSize', mkSZin)
        end
    end

% -------------
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
    
% ------------------------
% Plot diamonds at China Rock for referencing the other figure

%
ind_match_ref = find_ind_matchstring(instrtable.ChinaRock.locationID, ["B03", "B13"]);
    
%
plot(instrtable.ChinaRock.X(ind_match_ref), instrtable.ChinaRock.Y(ind_match_ref), 'kd', 'MarkerSize', 16, 'MarkerFaceColor', 'k')
%
plot(instrtable.ChinaRock.X(ind_match_ref(1)), instrtable.ChinaRock.Y(ind_match_ref(1)), 'yd', 'MarkerSize', 13, 'MarkerFaceColor', 'y')
plot(instrtable.ChinaRock.X(ind_match_ref(2)), instrtable.ChinaRock.Y(ind_match_ref(2)), 'rd', 'MarkerSize', 13, 'MarkerFaceColor', 'r')
    
% ------------------------
% Plot arrow to denote where photo was taken from
%
plt_arrow(haxs_CHR_land, 1.2, 0, 140, -16, -30)


% --------------------------------------------
% Label instrument sites

%
lblsiteFS = 18;
%
text(haxs_CHR_land, -795, 10, 'B03', 'Color', 'y', 'FontSize', lblsiteFS)
text(haxs_CHR_land, -200, -80, 'B13', 'Color', 'r', 'FontSize', lblsiteFS)

%
text(haxs_ASL_land, -845, -65, 'X01', 'Color', 'y', 'FontSize', lblsiteFS)


% --------------------------------------------
% Legend on top of China Rock panel

%
fill(haxs_CHR_land, [-800, -800, -450, -450, -800]-60, [-450, -220, -220, -450, -450], 'w')


%
xdotleg = -780 -60;
ydotleg = -245 - 60.*(0:3);
%
for i = 1:length(ydotleg)
    plot(haxs_CHR_land, xdotleg, ydotleg(i), '.k', 'MarkerSize', 36)
end
plot(haxs_CHR_land, xdotleg, ydotleg(1), '.b', 'MarkerSize', 28)
plot(haxs_CHR_land, xdotleg, ydotleg(2), '.r', 'MarkerSize', 28)
plot(haxs_CHR_land, xdotleg, ydotleg(3), '.y', 'MarkerSize', 28)
plot(haxs_CHR_land, xdotleg, ydotleg(4), '.g', 'MarkerSize', 28)
%
txtleg = 10;
%
text(haxs_CHR_land, xdotleg + 20, ydotleg(1), 'pressure', 'FontSize', txtleg)
text(haxs_CHR_land, xdotleg + 20, ydotleg(2), 'ADCP', 'FontSize', txtleg)
text(haxs_CHR_land, xdotleg + 20, ydotleg(3), 'Spotter', 'FontSize', txtleg)
text(haxs_CHR_land, xdotleg + 20, ydotleg(4), 'Spotter + pressure', 'FontSize', txtleg)

% -------
%
haxs_ASL_land.Visible = 'off';
haxs_CHR_land.Visible = 'off';


% --------------------------------------------
% Letter labeling

%
ltlblFS = 22;
%
text(haxs_a, 3310, 1700, 'a)', 'FontSize', ltlblFS)
text(haxs_c, 10, 180, 'b)', 'FontSize', ltlblFS)
%
text(haxs_d, -855, 65, 'c)', 'FontSize', ltlblFS)
text(haxs_e, -855, 410, 'd)', 'FontSize', ltlblFS)

%
text(haxs_d, -785, 65, 'Asilomar', 'FontSize', ltlblFS)
text(haxs_e, -780, 410, 'China Rock', 'FontSize', ltlblFS)



%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure01.pdf'), 'Resolution', 300)
% % %
% % exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
% %                               'figure01.png'), 'Resolution', 300)
    