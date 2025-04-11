%% Make Figure 9 -- plot dF/dx vs. dFx/dx

clear
close all


%% Load L4 data

%
load(fullfile(paper_directory(), 'data', 'level_4', 'roxsi_data_L4.mat'))


%% Define pairs for fe calculation

% Pairs of ADCPs
list_pairs = ["B11", "B13"; ...
              "B13", "B15"];


%% Run pairwise fe calculations

%
tic
disp('--- Doing pairwise fe estimate ---')

%
[dataL4, pwdissipation] = pairwise_fe_getdatapair(dataL4, list_pairs);

%
pwdissipation = pairwise_fe_calculate(pwdissipation);


%% Get standard deviation for each pair

pwdissipation = pairwise_fe_get_zmsl_stddev(pwdissipation);


%% Now apply QC and redo analysis

%
pwdissipation_flagged = pairwise_fe_QCflag(pwdissipation);
%
pwdissipation_flagged = pairwise_fe_calculate(pwdissipation_flagged);

%
pwfeQCtable_flagged = pairwise_fe_createQCtable(pwdissipation_flagged);


%% Overwrite pwdissipation with data structure after QC

pwdissipation = pwdissipation_flagged;


%% Compute regresssion

%
for i = 1:length(pwdissipation.data)

    %
    x_var = -pwdissipation.data(i).fromEcg.dFluxdx;
    y_var = -pwdissipation.data(i).fromFx.dFluxdx;
    
    %
    lok = ~isnan(x_var) & ~isnan(y_var);
    
    %
    x_var_ok = x_var(lok);
    y_var_ok = y_var(lok);
    %
    x_var_ok = x_var_ok(:);
    y_var_ok = y_var_ok(:);
    
    %
    mcoefs_aux = regress(y_var_ok, x_var_ok);
    %
    rcoef_aux = corr(x_var_ok, y_var_ok);
    
    %
    pwdissipation.data(i).EcgvsFx.rcoef = rcoef_aux;
    pwdissipation.data(i).EcgvsFx.linecoefs = mcoefs_aux;

end



%% Plot figure

%
xy_axs_A = [-8, 120];
xy_axs_B = [-3, 60];

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.39, 0.1483, 0.18, 0.4908];
%
haxs = makeSubPlots(0.08, 0.04, 0.1, ...
                    0.04, 0.15, 0.075, 1, 2);
hold(haxs, 'on')

    %
    mkSZ = 18;

    %
    plot(haxs(1), -pwdissipation.data(1).fromEcg.dFluxdx, ...
                  -pwdissipation.data(1).fromFx.dFluxdx, '.k', 'MarkerSize', mkSZ)
    %
    plot(haxs(2), -pwdissipation.data(2).fromEcg.dFluxdx, ...
                  -pwdissipation.data(2).fromFx.dFluxdx, '.k', 'MarkerSize', mkSZ)

	%
    plot(haxs(1), xy_axs_A, xy_axs_A, '--k', 'LineWidth', 1.5)
    plot(haxs(2), xy_axs_B, xy_axs_B, '--k', 'LineWidth', 1.5)
    
    %
    line_A = pwdissipation.data(1).EcgvsFx.linecoefs;
    line_B = pwdissipation.data(2).EcgvsFx.linecoefs;

    % Crossing through the origin
    plot(haxs(1), xy_axs_A, line_A.*xy_axs_A, '-r', 'LineWidth', 1.5)
    plot(haxs(2), xy_axs_B, line_B.*xy_axs_B, '-r', 'LineWidth', 1.5)                    
                        
%
set(haxs, 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
set(haxs, 'DataAspectRatio', [1, 1, 1])
set(haxs(1), 'XLim', xy_axs_A, 'YLim', xy_axs_A)
set(haxs(2), 'XLim', xy_axs_B, 'YLim', xy_axs_B)
%
set(haxs(1), 'XTick', 0:25:125, 'YTick', 0:25:125)
set(haxs(2), 'XTick', 0:20:60, 'YTick', 0:20:60)
%
set(haxs(1), 'XTickLabelRotation', 0)

%
site_1_A = char(pwdissipation.locationID(1, 1));
site_1_B = char(pwdissipation.locationID(1, 2));
%
site_2_A = char(pwdissipation.locationID(2, 1));
site_2_B = char(pwdissipation.locationID(2, 2));

%
text(haxs(1), -6, 108, ['a) ' site_1_A '-' site_1_B], 'Interpreter', 'Latex', 'FontSize', 22)
text(haxs(2), -2, 54, ['b) ' site_2_A '-' site_2_B], 'Interpreter', 'Latex', 'FontSize', 22)

%
lblFS = 16;
%
xlabel(haxs(end), '$-\mathrm{d}F/\mathrm{d} x$ [W m$^{-2}$]', 'Interpreter', 'Latex', 'FontSize', lblFS)
ylabel(haxs(1), '$-\mathrm{d}F_x/\mathrm{d} x$  [W m$^{-2}$]', 'Interpreter', 'Latex', 'FontSize', lblFS)
ylabel(haxs(2), '$-\mathrm{d}F_x/\mathrm{d} x$  [W m$^{-2}$]', 'Interpreter', 'Latex', 'FontSize', lblFS)


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure09.pdf'), 'Resolution', 300)

