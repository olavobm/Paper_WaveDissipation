%% Make figure 3 -- overview plot with timeseries of wave conditions.

clear
close all


%% Load L4 data

%
dataL4 = load(fullfile(paper_directory(), 'data', 'level_4', 'roxsi_data_L4.mat'));
dataL4 = dataL4.dataL4;


%% Choose which locations to plot

%
pltID = ["B03", "B13"];

%
indID = find_ind_matchstring(dataL4.locationID, pltID);


%% Time limits for the plot

%
dtime_lims = [datetime(2022, 06, 16, 00, 00, 00), ...
              datetime(2022, 07, 20, 00, 00, 00)];
dtime_lims.TimeZone = 'America/Los_Angeles';

                          
%% Make figure

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.15, 0.0967, 0.25, 0.5];
%
haxs = makeSubPlots(0.15, 0.08, 0.1, ...
                    0.025, 0.1, 0.02, 1, 3);
hold(haxs, 'on')

    %
    plot(haxs(1), dataL4.dtime, dataL4.seaswell.Hs(:, indID), 'LineWidth', 2)
    plot(haxs(2), dataL4.dtime, 1./dataL4.seaswell.meanfrequency(:, indID), 'LineWidth', 2)
    plot(haxs(3), dataL4.dtime, dataL4.seaswell.meandir1(:, indID), 'LineWidth', 2)
    
    %
    plot(haxs(3), dtime_lims, [0, 0], '-k')
    
    
    
% ---------------------------------------
%
hleg = legend(haxs(2), [char(pltID(1)) ', $\overline{h} = ' num2str(mean(dataL4.bottomdepth(:, indID(1)), 1, 'omitnan'), '%.1f') '~\mathrm{m}$'], ...
                       [char(pltID(2)) ', $\overline{h} = ' num2str(mean(dataL4.bottomdepth(:, indID(2)), 1, 'omitnan'), '%.1f') '~\mathrm{m}$'], ...
                       'Interpreter', 'Latex', 'Location', 'SouthWest');
    hleg.FontSize = 16;
    hleg.Position = [0.555, 0.5779, 0.3643, 0.0978];
    
    
% ---------------------------------------
%
set(haxs, 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
%
set(haxs(1), 'YLim', [0, 2])
set(haxs(2), 'YLim', [5.5, 13])
set(haxs(3), 'YLim', [-40, 40])
%
set(haxs(1), 'YTick', [0, 0.5, 1, 1.5, 2])
%
set(haxs, 'XLim', dtime_lims)
set(haxs(1:end-1), 'XTickLabel', [])

% ---------------------------------------
%
lblFS = 18;
%
ylabel(haxs(1), '$H_s$ [m]', 'Interpreter', 'Latex', 'FontSize', lblFS)
ylabel(haxs(2), '$T_{\mathrm{mean}}$ [s]', 'Interpreter', 'Latex', 'FontSize', lblFS)
ylabel(haxs(3), '$\theta_{\mathrm{mean}}$ [degrees]', 'Interpreter', 'Latex', 'FontSize', lblFS)


% ---------------------------------------
% Different letter labeling

%
xtxtplt_1 = dataL4.dtime(18);
%
txtFS = 30;
txtFS_letter = txtFS - 8;
%
text(haxs(1), xtxtplt_1, 1.75, 'a)', 'Interpreter', 'Latex',  'FontSize', txtFS_letter)
text(haxs(2), xtxtplt_1, 12.1, 'b)', 'Interpreter', 'Latex',  'FontSize', txtFS_letter)
text(haxs(3), xtxtplt_1, 30, 'c)', 'Interpreter', 'Latex',  'FontSize', txtFS_letter)


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure03.pdf'), 'Resolution', 300)

                          