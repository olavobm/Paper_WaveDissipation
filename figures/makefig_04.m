%% Make Figure 4 -- cross-shore statistics

clear
close all


%% Load L4 data

%
dataAll = load(fullfile(paper_directory(), 'data', 'level_4', 'roxsi_data_L4.mat'));
dataAll = dataAll.dataL4;


%% Select instrument sites to plot

%
pick_sites = ["B03"; "B05"; "B06"; "B07"; ...
              "B09"; "B10"; "B11"; "B12"; "B13"; ...
              "B14"; "B15"; "B16"; ...
              "A01"; "A02"; "A04"; "A05"; "A06"; "A07"; ...
              "C01"; "C02"; ....
              "D01"; "D02"; ...
              "X01"; "X03"; "X04"; "X05"; "X06"; ...
              "X07"; "X08"; "X09"; "X10"; "X11"; ...
              "E01"; "E02"; ...
              "E03"; ...
              "E04"; "E05"; "E06"; "E08"; "E11"; "E12"; "E13"];

%
ind_match = find_ind_matchstring(dataAll.locationID, pick_sites);

%
Nsites = length(pick_sites);

          
%% Check common time of data available

%
lok_common = ~isnan(mean(dataAll.bottomdepth(:, ind_match), 2)) & ...
             ~isnan(mean(dataAll.seaswell.Hs(:, ind_match), 2));

%
hfig = figure;
    plot(dataAll.dtime, lok_common)
    ylim([-0.1, 1.1])
%
title('data availability (1 = ok; 0 = NaN)')

close(hfig);


%% Calculate time-mean and quantile ranges
  
%
h_avg = mean(dataAll.bottomdepth(lok_common, ind_match), 1);
Hs_avg = mean(squeeze(dataAll.seaswell.Hs(lok_common, ind_match)), 1);
Ecg_avg = mean(squeeze(dataAll.seaswell.Ecg(lok_common, ind_match)), 1);
%
Hs_q1 = quantile(squeeze(dataAll.seaswell.Hs(lok_common, ind_match)), 0.25, 1);
Hs_q2 = quantile(squeeze(dataAll.seaswell.Hs(lok_common, ind_match)), 0.75, 1);
%
Ecg_q1 = quantile(squeeze(dataAll.seaswell.Ecg(lok_common, ind_match)), 0.25, 1);
Ecg_q2 = quantile(squeeze(dataAll.seaswell.Ecg(lok_common, ind_match)), 0.75, 1);


%% Compute mean and quantiles of normalized Flux
% (flux at China Rock and Asilomar are
% normalized by the corresponding offshore instruments)

%
allEcg = squeeze(dataAll.seaswell.Ecg(:, ind_match));
allEcg_norm = allEcg;
%
allEcg_norm(~lok_common, :) = NaN;     % make all timeseries consistent

%
ind_ASL = find(strcmp(dataAll.site(ind_match), "Asilomar"));
ind_CHR = find(strcmp(dataAll.site(ind_match), "ChinaRock"));

% Normalize by the corresponding offshore measurements
allEcg_norm(:, ind_ASL) = allEcg_norm(:, ind_ASL) ./  allEcg_norm(:, ind_ASL(1));
allEcg_norm(:, ind_CHR) = allEcg_norm(:, ind_CHR) ./  allEcg_norm(:, ind_CHR(1));

%
allEcg_norm_avg = mean(allEcg_norm(lok_common, :), 1);
allEcg_norm_q1 = quantile(allEcg_norm(lok_common, :), 0.25, 1);
allEcg_norm_q2 = quantile(allEcg_norm(lok_common, :), 0.75, 1);


%% Run TG83 without bottom friction

%
dx = 0.5;
xvec = -1000:dx:-5;
%
h = -tand(atand(1/40)).*xvec;

%
freq0 = 1/7.8799;
%
Hrms0 = 0.9565/sqrt(2);
%
angle0 = 0;
%
gamma = 0.42;
B = 1;

%
modeloutput = wave_xshoretransformation(xvec, h, dx, ...
                                        freq0, Hrms0, angle0, ...
                                        gamma, B);

%% Make figure

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.36, 0.1, 0.275, 0.55];
%
haxs = makeSubPlots(0.15, 0.05, 0.1, ...
                    0.08, 0.1, 0.02, 1, 2);
hold(haxs, 'on')

% ----------------------------------------------------
% Plot symbols as function of bottom depth

%
mkSZ = 32;

    %
    plot(haxs(1), h_avg, Hs_avg, '.k', 'MarkerSize', mkSZ)
    %
    plot(haxs(2), h_avg, allEcg_norm_avg, '.k', 'MarkerSize', mkSZ)
    
    
    %
    for i = 1:Nsites
        % ---------------------------------
        %
        dHs_1_aux = Hs_avg(i) - Hs_q1(i);
        dHs_2_aux = Hs_q2(i) - Hs_avg(i);
        %
        plot(haxs(1), h_avg(i).*[1, 1], [Hs_q1(i), Hs_q2(i)], '-k')
        
        %
        plot(haxs(2), h_avg(i).*[1, 1], [allEcg_norm_q1(i), allEcg_norm_q2(i)], '-k')
    end


% ----------------------------------------------------
% Plot TG83 model

hp_TG83_Hs = plot(haxs(1), modeloutput.h, sqrt(2)*modeloutput.Hrms, '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2);
plot(haxs(2), modeloutput.h, modeloutput.Ecgx./modeloutput.Ecgx(1), '-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2)

    
% ----------------------------------------------------
% Plot B03 and B13 in different symbols

%
plthlID = ["B03", "B13"];
%
ID_allmatch = dataAll.locationID(ind_match);
instr_allmatch = dataAll.instrument(ind_match);

%
for i = 1:length(plthlID)
    %
    indpltsymbol = find(strcmp(ID_allmatch, plthlID(i)));
    %
    if strcmp(instr_allmatch(indpltsymbol), "Spotter")
        clrplt = 'y';
    elseif strcmp(instr_allmatch(indpltsymbol), "SoloD")
        clrplt = 'b';
    elseif strcmp(instr_allmatch(indpltsymbol), "Signature") || ...
           strcmp(instr_allmatch(indpltsymbol), "Aquadopp") 
        clrplt = 'r';
    end

    %
    plot(haxs(1), h_avg(indpltsymbol), Hs_avg(indpltsymbol), 'd', ...
                  'MarkerSize', 16, 'Color', 'k', 'MarkerFaceColor', 'k');
    plot(haxs(1), h_avg(indpltsymbol), Hs_avg(indpltsymbol), 'd', ...
                  'MarkerSize', 12, 'Color', clrplt, 'MarkerFaceColor', clrplt);
    %
    plot(haxs(2), h_avg(indpltsymbol), allEcg_norm_avg(indpltsymbol), 'd', ...
                  'MarkerSize', 16, 'Color', 'k', 'MarkerFaceColor', 'k');
    plot(haxs(2), h_avg(indpltsymbol), allEcg_norm_avg(indpltsymbol), 'd', ...
                  'MarkerSize', 12, 'Color', clrplt, 'MarkerFaceColor', clrplt);       
              
end
   

% ----------------------------------------------------
% Set properties of axes

%
set(haxs, 'FontSize', 16, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
set(haxs, 'XDir', 'reverse', 'XLim', [0, 24])
%
set(haxs(1), 'YLim', [0, 1.3])
set(haxs(2), 'YLim', [0, 1.4])

%
set(haxs, 'XTick', 0:2:30)
set(haxs(1:end-1), 'XTickLabel', [])



% ----------------------------------------------------
% Create legend manually

%
hrect = fill(haxs(1), [23.5, 23.5, 9.8, 9.8, 23.5], [0.04, 0.4, 0.4, 0.04, 0.04], 'w');

%
yposleg = 0.33 - 0.11.*(0:2);

%
plot(haxs(1), 22.6, yposleg(1), 'dk', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
plot(haxs(1), 22.6, yposleg(2), 'dk', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
%
plot(haxs(1), 22.6, yposleg(1), 'dy', 'MarkerFaceColor', 'y', 'MarkerSize', 8)
plot(haxs(1), 22.6, yposleg(2), 'dr', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
%
plot(haxs(1), [23.2, 22], yposleg(3).*[1, 1], 'Color', hp_TG83_Hs.Color, 'LineWidth', 3)

%
text(haxs(1), 21.75, yposleg(1), 'B03', 'FontSize', 14)
text(haxs(1), 21.75, yposleg(2), 'B13', 'FontSize', 14)
text(haxs(1), 21.75, yposleg(3), '1D model (no bottom friction)', 'FontSize', 14)

% ----------------------------------------------------
% Add x/y labels

%
xlabel(haxs(end), '$\overline{h}$ [m]', 'Interpreter', 'Latex', 'FontSize', 22)
%
ylabel(haxs(1), '$H_s$ [m]', 'Interpreter', 'Latex', 'FontSize', 22)
ylabel(haxs(2), '$F / F_0$', 'Interpreter', 'Latex', 'FontSize', 22)


% ----------------------------------------------------
% Letter labelling

%
xlbltxt = 23.75;
lblFS = 26;
%
text(haxs(1), xlbltxt, 1.175, 'a)', 'Interpreter', 'Latex', 'FontSize', lblFS)
text(haxs(2), xlbltxt, 1.27, 'b)', 'Interpreter', 'Latex', 'FontSize', lblFS)


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure04.pdf'), ...
                     'Resolution', 300)
