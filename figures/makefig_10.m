%% Make Figure 10 -- fe vs. Ab/sigma_h (ROXSI and literature)

clear
close all


%%
% -----------------------------------------------
% ------------------ LOAD DATA ------------------
% -----------------------------------------------

%% Load bin-averaged results

%
data_file = fullfile(paper_directory(), 'data', 'level_5', 'bin_averaged_data.mat');
%
feROXSI2022 = load(data_file);
feROXSI2022 = feROXSI2022.feROXSI2022;


%%
% -----------------------------------------------
% --------- DEFINE A/k GRID AND COMPUTE ---------
% --------------- ROXSI2022 CURVE ---------------
% -----------------------------------------------


%% Set-up A/k grid

%
Aonsigmah_grid = logspace(log10(0.05), log10(200), 100);
kNonsigmah = 4;
AonkN_grid = Aonsigmah_grid./kNonsigmah;

%
allresults.Aonsigmah = Aonsigmah_grid;
allresults.AonkN = AonkN_grid;


%% Power-law from ROXSI 2022

% Regression from the bin means
x_aux = log10(feROXSI2022.Abonsigma_grid);
y_aux = log10(feROXSI2022.fe_binmean_sigmahvar);
%
mcoefs_aux = regress(y_aux, [ones(length(x_aux), 1), x_aux]);

%
allresults.ROXSI22.powerlaw_coefs = [(10.^mcoefs_aux(1)), mcoefs_aux(2)];

%
allresults.ROXSI22.fe = (10.^mcoefs_aux(1)) .* (Aonsigmah_grid.^mcoefs_aux(2));


%%
% -----------------------------------------------
% ---------- COMPUTE CURVES OF fe FROM ----------
% -------------- PARAMETERIZATIONS --------------
% -----------------------------------------------


%% Jonsson and Carlsen (1976) equation

%
JC76_RHS = 0.2 + log10(AonkN_grid);

%
JC76_findLHS = @(f, RHS) 1./(4*sqrt(f)) + log10(1./(4*sqrt(f))) - RHS;
%
fw_JC76 = NaN(size(JC76_RHS));

%
for i = 1:length(JC76_RHS)
    %
    fw_JC76(i) = fzero(@(f) JC76_findLHS(f, JC76_RHS(i)), [1e-5, 5000]);
end

%
allresults.JC76.fw = fw_JC76;
allresults.JC76.fe = allresults.JC76.fw;


%% Madsen (1994) -- approximation equation 32, with coefficients
% for large roughness, supposedly valid for A/kN >= 0.2

%
a1_M94_cte = 7.02;
a2_M94_cte = -0.078;
a3_M94_cte = -8.82;

%
fw_M94 = @(Aonk) exp(a1_M94_cte .* (Aonk).^a2_M94_cte + a3_M94_cte);

% Angle in degrees
theta_M94 = @(Aonk) 33 - 6 * log10(Aonk);
%
allresults.costheta_M4 = cosd(theta_M94(AonkN_grid));

%
allresults.M94.fw = fw_M94(AonkN_grid);
allresults.M94.fe = allresults.M94.fw.*allresults.costheta_M4;


%% Rogers et al. (2016) parameterization

%
R16_f = NaN(length(AonkN_grid), 1);
%
l_bigTH = (AonkN_grid > 0.0369);
R16_f(l_bigTH) = exp(5.213 * (AonkN_grid(l_bigTH)).^(-0.194) + -5.977);
%
R16_f(~l_bigTH) = 50;

%
allresults.R16.fw = R16_f;
allresults.R16.fe = allresults.R16.fw;


%%
% -----------------------------------------------
% ---------- DATA FROM THE LITERATURE -----------
% -----------------------------------------------


%% Lowe et al. (2005)

%
Ab_Lowe = [0.3; 0.27; 0.23]./[1.02; 1.01; 0.98];
%
allresults.L05.Aonsigmah = Ab_Lowe./0.036;
allresults.L05.fe = [0.23; 0.24; 0.24];


%% Lentz et al. (2016)

%
allresults.L16.Aonsigmah = [0.03; 0.05; 0.07; 0.09; 0.11; 0.13; 0.16; 0.215; 0.3];
allresults.L16.Aonsigmah = allresults.L16.Aonsigmah./0.13;
%
allresults.L16.fe = [4.2338; 2.9213; 2.3959; 1.959; 1.5956; 1.4613; 1.3337; 1.1875; 0.9058];


%% Gon et al. (2020)

%
allresults.G20.Aonsigmah = [0.0983; 0.1736; 0.2651; 0.3869; 0.4915; 0.6084; 0.7074];
allresults.G20.fe = [31.56; 15.21; 10.41; 5.13; 4.11; 2.50; 2.06];

%
x_aux = log10(allresults.G20.Aonsigmah(:));
y_aux = log10(allresults.G20.fe(:));
%
mcoefs_aux = regress(y_aux, [ones(length(x_aux), 1), x_aux]);


%% Sous et al. (2023)

%
data_Sous.Zone34.Abonsigma_grid = [2; 3; 4; 5; 6; 7; 8];
data_Sous.Zone34.fe_binmean = [0.774248170392414; ...
                               0.585161629822249; ...
                               0.482443036347996; ...
                               0.432033897855224; ...
                               0.399898539824951; ...
                               0.391455993979168; ...
                               0.374442400805494];
%
data_Sous.Zone45.Abonsigma_grid = [1.5; 2; 2.5; 3.5];
data_Sous.Zone45.fe_binmean = [1.258673930054496; ...
                               1.034097527988925; ...
                               0.887875368349446; ...
                               0.735483918948416];
%
data_Sous.Zone56.Abonsigma_grid = [0.4; 0.6; 0.8];
data_Sous.Zone56.fe_binmean = [5.429890861132564; ...
                               4.494433556698649; ...
                               3.534580800451720];

%
allresults.S23.R1.Aonsigmah = data_Sous.Zone34.Abonsigma_grid(:);
allresults.S23.R1.fe = data_Sous.Zone34.fe_binmean(:);
%
allresults.S23.R2.Aonsigmah = data_Sous.Zone45.Abonsigma_grid(:);
allresults.S23.R2.fe = data_Sous.Zone45.fe_binmean(:);
%
allresults.S23.R3.Aonsigmah = data_Sous.Zone56.Abonsigma_grid(:);
allresults.S23.R3.fe = data_Sous.Zone56.fe_binmean(:);


%%
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% --------------------------- MAKE FIGURE ---------------------------
% -------------------------------------------------------------------
% -------------------------------------------------------------------


%% Make figure

%
lnWd = 1.5;

% -----------------------------------------------
%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.5    0.2492    0.28    0.5025];
%
haxs = makeSubPlots(0.18, 0.1, 0.1, ...
                    0.1, 0.15, 0.1, 1, 1);
hold(haxs, 'on')
 
    % -----------------------------------------------
    % Plot theory/semi-empirical relationships

    %
    hp_JC76 = plot(allresults.Aonsigmah, 0.75*allresults.JC76.fe,  'LineWidth', lnWd);
    hp_M94 = plot(allresults.Aonsigmah, 0.875*allresults.M94.fe, 'LineWidth', lnWd);

    %
    hp_R16curve = plot(allresults.Aonsigmah, 0.75*allresults.R16.fe, 'LineWidth', lnWd);

    
    % -----------------------------------------------
    % Plot data from the literature

    %
    mkSz = 38;
    %
    hp_L05 = plot(allresults.L05.Aonsigmah, (0.7/0.8)*allresults.L05.fe, '.', 'MarkerSize', mkSz);
    %
    hp_L16 = plot(allresults.L16.Aonsigmah, 0.5*allresults.L16.fe, '.', 'Color', [0.8500    0.3250    0.0980], 'MarkerSize', mkSz);
    %
    hp_G20 = plot(allresults.G20.Aonsigmah, 0.5*allresults.G20.fe, '.', 'MarkerSize', mkSz);
    
    %
    hp_S23 = plot(allresults.S23.R1.Aonsigmah, (0.7/0.8)*allresults.S23.R1.fe, '.', 'MarkerSize', mkSz);
    plot(allresults.S23.R2.Aonsigmah, (0.7/0.8)*allresults.S23.R2.fe, '.', 'Color', hp_S23.Color, 'MarkerSize', mkSz);
    plot(allresults.S23.R3.Aonsigmah, (0.7/0.8)*allresults.S23.R3.fe, '.', 'Color', hp_S23.Color, 'MarkerSize', mkSz);
   

    % -----------------------------------------------
    % Plot ROXSI data
    
    %
    for i = 1:length(feROXSI2022.Abonsigma_grid)
        
        %
        plot(feROXSI2022.Abonsigma_grid(i).*[1, 1], ...
             feROXSI2022.fe_quant_sigmahvar(i, :), '-k', 'LineWidth', 2)
        
    end
    
    %
    hp_ROXSI2022data = plot(feROXSI2022.Abonsigma_grid, ...
                            feROXSI2022.fe_binmean_sigmahvar, '.k', 'MarkerSize', 48);
    
     
	% -----------------------------------------------
    % Plot ROXSI power law
    %
    hp_ROXSIscaling = plot(allresults.Aonsigmah, allresults.ROXSI22.fe, '-k', 'LineWidth', 2.5);
     
    

% -----------------------------------------------
% Create legend manually

%
width_factor = 4.5;
%
xlA = 0.6;
xlB = 3;

%
hrect_1 = fill([xlA, xlA, xlA*width_factor, xlA*width_factor, xlA], [5, 55, 55, 5, 5], 'w');
%
hrect_2 = fill([xlB, xlB, xlB*width_factor, xlB*width_factor, xlB], [10, 55, 55, 10, 10], 'w');

%
y_pos_aux = 40;
y_pos_aux = y_pos_aux .* 0.55.^(0:1:3);
%
x_pos_aux = xlA.*[1.05, 1.5];

lnWdleg = 4;
%
plot(x_pos_aux, y_pos_aux(1).*[1, 1], '-', 'LineWidth', lnWdleg, 'Color', hp_JC76.Color);
plot(x_pos_aux, y_pos_aux(2).*[1, 1], '-', 'LineWidth', lnWdleg, 'Color', hp_M94.Color);
plot(x_pos_aux, y_pos_aux(3).*[1, 1], '-', 'LineWidth', lnWdleg, 'Color', hp_R16curve.Color);
plot(x_pos_aux, y_pos_aux(4).*[1, 1], '-', 'LineWidth', lnWdleg, 'Color', hp_ROXSIscaling.Color);

%
text(x_pos_aux(2).*1.05, y_pos_aux(1), 'JC76', 'Interpreter', 'Latex', 'FontSize', 16);
text(x_pos_aux(2).*1.05, y_pos_aux(2), 'M94', 'Interpreter', 'Latex', 'FontSize', 16);
text(x_pos_aux(2).*1.05, y_pos_aux(3), 'R16', 'Interpreter', 'Latex', 'FontSize', 16);
text(x_pos_aux(2).*1.05, y_pos_aux(4), 'ROXSI22', 'Interpreter', 'Latex', 'FontSize', 16);

%
y_pos_aux = 40;
y_pos_aux = y_pos_aux .* 0.55.^(0:1:3);

%
mkSZleg = 42;

%
x_pos_aux = xlB.*1.15;

%
plot(x_pos_aux, y_pos_aux(1), '.', 'Color', hp_L05.Color, 'MarkerSize', mkSZleg);
plot(x_pos_aux.*2.2, y_pos_aux(1), '.', 'Color', hp_L16.Color, 'MarkerSize', mkSZleg);
plot(x_pos_aux.*2.2, y_pos_aux(2), '.', 'Color', hp_G20.Color, 'MarkerSize', mkSZleg);
plot(x_pos_aux, y_pos_aux(2), '.', 'Color', hp_S23.Color, 'MarkerSize', mkSZleg);
plot(x_pos_aux, y_pos_aux(3), '.', 'Color', hp_ROXSI2022data.Color, 'MarkerSize', mkSZleg);

%
txtSZleg = 16;
%
text(x_pos_aux.*1.2, y_pos_aux(1), 'L05', 'Interpreter', 'Latex', 'FontSize', txtSZleg);
text(x_pos_aux.*2.5, y_pos_aux(1), 'L16', 'Interpreter', 'Latex', 'FontSize', txtSZleg);
text(x_pos_aux.*2.5, y_pos_aux(2), 'G20', 'Interpreter', 'Latex', 'FontSize', txtSZleg);
text(x_pos_aux.*1.2, y_pos_aux(2), 'S23', 'Interpreter', 'Latex', 'FontSize', txtSZleg);
text(x_pos_aux.*1.2, y_pos_aux(3), 'ROXSI22', 'Interpreter', 'Latex', 'FontSize', txtSZleg);


% -----------------------------------------------
% Set axes

%
set(gca, 'FontSize', 16, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
%
set(gca, 'XScale', 'log', 'YScale', 'log')
set(gca, 'XLim', [0.09, 15], 'YLim', [1e-1, 60])
set(gca, 'XTick', [0.1, 1, 10])

%
xlabel('$A_b/\sigma_h$', 'Interpreter', 'Latex', 'FontSize', 28)
ylabel('$f_e$', 'Interpreter', 'Latex', 'FontSize', 28)


%% Save figure

%
dir_output = fullfile(paper_directory(), 'figures');

%
exportgraphics(hfig, fullfile(dir_output, 'figure10.pdf'), 'Resolution', 300)


