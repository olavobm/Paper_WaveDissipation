%% Georectify Monterey Peninsula image


clear
close all


%% Load image

%
photo_MP_georect = imread(fullfile(paper_directory(), 'figures', 'MP_googleearth_highres.jpg'));


%% 3 reference locations

% (at NPS)
reflocs.NPS.x = 4818.442017271023;
reflocs.NPS.y = 2558.835321551710;
%
reflocs.NPS.latitude = 36.594508333333;
reflocs.NPS.longitude = -121.87436111111;
%     36 35 40.23
%     -121 -52 - -27.70

% near Point Pinos
reflocs.ptPinos.x = 4032.278812326220;
reflocs.ptPinos.y = 1891.874027034726;
%
reflocs.ptPinos.latitude = 36.636505555555556;
reflocs.ptPinos.longitude = -121.934652777777778;
%     36 38 11.42
%     -121 -56 - 4.75


% near Sunset Point (Pebble Beach)
reflocs.sunsetpoint.x = 3618.910752646380 ;
reflocs.sunsetpoint.y = 2970.383339165542;
%
reflocs.sunsetpoint.latitude = 36.569711111;
reflocs.sunsetpoint.longitude = -121.967622222222222;
%     36 34 10.96
%     -121 -58 - 3.44


%%
% ----------------------------------------------
% ---- CONVERSION FROM (X, Y) TO (LAT, LON) ----
% ----------------------------------------------


%% Axis limits
    
%
xaxislims = [3300, 5000];
yaxislims = [1600, 3500];


%% Uses axis limits to get relative coordinates for leat-squares fit

%
x0 = mean(xaxislims);
y0 = mean(yaxislims);

%
xprime = [reflocs.NPS.x; reflocs.ptPinos.x; reflocs.sunsetpoint.x] - x0;
yprime = [reflocs.NPS.y; reflocs.ptPinos.y; reflocs.sunsetpoint.y] - y0;

%
latref = [reflocs.NPS.latitude; reflocs.ptPinos.latitude; reflocs.sunsetpoint.latitude];
lonref = [reflocs.NPS.longitude; reflocs.ptPinos.longitude; reflocs.sunsetpoint.longitude];


%%

%
coefsgeorect.x0 = x0;
coefsgeorect.y0 = y0;

%
coefsgeorect.data = reflocs;


%% Do least-squares for longitude

%
Gm = [ones(size(xprime)), xprime];
%
Gm_aux = (Gm.' * Gm) \ Gm.';

%
coefsgeorect.xprime_to_lon = Gm_aux * lonref;


%% Do least-squares for latitude

%
Gm = [ones(size(yprime)), yprime];
%
Gm_aux = (Gm.' * Gm) \ Gm.';

%
coefsgeorect.yprime_to_lat = Gm_aux * latref;


%% Create equation handle


%
coefsgeorect.fcn_lon = @(x) coefsgeorect.xprime_to_lon(1) + coefsgeorect.xprime_to_lon(2).*(x - coefsgeorect.x0);
coefsgeorect.fcn_lat = @(y) coefsgeorect.yprime_to_lat(1) + coefsgeorect.yprime_to_lat(2).*(y - coefsgeorect.y0);


%%
% ----------------------------------------------
% ---- CONVERSION FROM (LAT, LON) TO (X, Y) ----
% ----------------------------------------------

%

%%

%
coefsgeorect.fcn_x = @(lon) ((lon - coefsgeorect.xprime_to_lon(1))./coefsgeorect.xprime_to_lon(2)) + coefsgeorect.x0;
coefsgeorect.fcn_y = @(lat) ((lat - coefsgeorect.yprime_to_lat(1))./coefsgeorect.yprime_to_lat(2)) + coefsgeorect.y0;



% %%
% return

%%
% --------------------------------------------------------------------
% --------------------------------------------------------------------
% ------------------- FIGURES TO CHECK THE RESULTS -------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------


%% Similar as above, but with Google Earth image on the left

%
lat_plt = 36 + (35./60);
lon_plt = -121 -(56./60);

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.3186    0.1842    0.3332    0.5358];
%
haxs_a = axes('Position', [0.1, 0.1, 0.8, 0.8]);   % Google Earth/MP


    % ----------------------------------------------------
    % Google Earth with MP
    
    %
    image(haxs_a, photo_MP_georect)
    %
    hold(haxs_a, 'on')
    %
    set(haxs_a, 'DataAspectRatio', [1, 1, 1])

    %
    plot(haxs_a, coefsgeorect.fcn_x(lon_plt), ...
                 coefsgeorect.fcn_y(lat_plt), '.r', 'MarkerSize', 42)
  
	%
    set(haxs_a, 'XLim', xaxislims, 'YLim', yaxislims)
 
    
%% Similar as above, but with longitude, latitude coordinates
    
%
xticks = 3400:200:5000;
yticks = 1600:200:3400;

% %
% lon_ticks = coefsgeorect.xprime_to_lon(1) + coefsgeorect.xprime_to_lon(2).*(xticks - coefsgeorect.x0);
% lat_ticks = coefsgeorect.yprime_to_lat(1) + coefsgeorect.yprime_to_lat(2).*(yticks - coefsgeorect.y0);


%
lon_ticks = coefsgeorect.fcn_lon(xticks);
lat_ticks = coefsgeorect.fcn_lat(yticks);


%
x_cell_label = cell(1, length(xticks));
for i = 1:length(xticks)
    x_cell_label{i} = num2str(lon_ticks(i), '%.2f');
end

y_cell_label = cell(1, length(yticks));
for i = 1:length(yticks)
    y_cell_label{i} = num2str(lat_ticks(i), '%.2f');
end

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.5    0.1842    0.3332    0.5358];
%
haxs_a = axes('Position', [0.1, 0.1, 0.8, 0.8]);   % Google Earth/MP


    % ----------------------------------------------------
    % Google Earth with MP
    
    %
    image(haxs_a, photo_MP_georect)
    %
    hold(haxs_a, 'on')
    %
    set(haxs_a, 'DataAspectRatio', [1, 1, 1])

%     %
%     plot(haxs_a, xyCHR(1), xyCHR(2), '+r', 'LineWidth', 4, 'MarkerFaceColor', 'r', 'MarkerSize', 22)
%     plot(haxs_a, xyASL(1), xyASL(2), '+r', 'LineWidth', 4, 'MarkerFaceColor', 'r', 'MarkerSize', 26)
    
    %
    set(haxs_a, 'XLim', xaxislims, 'YLim', yaxislims)
 
    %
    set(haxs_a, 'XTick', xticks, 'XTickLabel', x_cell_label, ...
                'YTick', yticks, 'YTickLabel', y_cell_label)

            
%% Save structure with georectification parameters

%
dir_output = mfilename('fullpath');
dir_output = fileparts(dir_output);

%
save(fullfile(dir_output, 'georect_MP.mat'), 'coefsgeorect')



            