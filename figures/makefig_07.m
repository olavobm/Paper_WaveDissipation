%% Make Figure 7 -- fe vs. Ab/sigma, with variable and constant sigma_h
%
% Also fit power-law and save in file feROXSI2022.mat.

clear
close all


%% Load L5 data

%
dir_data = fullfile(paper_directory(), 'data', 'level_5');

%
load(fullfile(dir_data, 'dataL5_withQC.mat'))


%% Keep only locations that passed QC

%
lpassed_QC = (pwfetable.pass_allQC);
%
ind_passed_QC = find(lpassed_QC);

%
pwfetable_passedQC = pwfetable(lpassed_QC, :);


%% Sort by correlation coefficient

%
[~, ind_sort] = sort(pwfetable_passedQC.rcoef, 'descend');

%
ind_passed_QC = ind_passed_QC(ind_sort);

%
pwfetable_passedQC = pwfetable_passedQC(ind_sort, :);


%% Get data and compute cumulative correlation
% of log_10(Ab/sigma_h) vs. log_10(fe)

%
N_aux = length(ind_passed_QC);

%
sigmah_avgsites = mean(pwfetable_passedQC.sigmah_avg);

%
datascaling.locID = strings(N_aux, 2);

%
datascaling.sigmah = NaN(N_aux, 1);
datascaling.sigmah_ref = sigmah_avgsites;

%
datascaling.ind_pairs = NaN(N_aux, 2);
datascaling.ind_pairs(1, 1) = 1;

%
datascaling.fe = [];
%
datascaling.Ab = [];
datascaling.Abonsigmah = [];
datascaling.Abonsigmah_ref = [];

%
datascaling.r_A_cumulative = NaN(N_aux, 1);
datascaling.r_B_cumulative = NaN(N_aux, 1);


%
for i = 1:length(ind_passed_QC)
    
    %% Get location ID
    
    datascaling.locID(i, 1) = pwdissipation.locationID(ind_passed_QC(i), 1);
    datascaling.locID(i, 2) = pwdissipation.locationID(ind_passed_QC(i), 2);
    
    %%
    
    %
    Ab_timeseries = pwdissipation.data(ind_passed_QC(i)).Ab_avg;
	%
    fe_timeseries = pwdissipation.data(ind_passed_QC(i)).fromEcg.fe_timeseries;
    %
    sigmah_aux = pwfetable_passedQC.sigmah_avg(i);
    
    %
    datascaling.sigmah(i) = sigmah_aux;
    %
    datascaling.fe = [datascaling.fe; fe_timeseries];
    datascaling.Ab = [datascaling.Ab; Ab_timeseries];
    
    %
    datascaling.Abonsigmah = [datascaling.Abonsigmah; (Ab_timeseries./sigmah_aux)];
    datascaling.Abonsigmah_ref = [datascaling.Abonsigmah_ref; (Ab_timeseries./sigmah_avgsites)];
    
    %
    datascaling.ind_pairs(i, 2) = length(datascaling.fe);
    %
    if i<length(ind_passed_QC)
        %
        datascaling.ind_pairs(i+1, 1) = datascaling.ind_pairs(i, 2) + 1;
    end
    
    
    %% Calculate correlation
    
    %
    x_A_var_aux = log10(datascaling.Abonsigmah);
    x_B_var_aux = log10(datascaling.Abonsigmah_ref);
    %
    y_var_aux = log10(datascaling.fe);
    
    %
    datascaling.r_A_cumulative(i) = corr(x_A_var_aux, y_var_aux, 'rows', 'complete');
    datascaling.r_B_cumulative(i) = corr(x_B_var_aux, y_var_aux, 'rows', 'complete');
    
    
end


%% Set pair number threshold and get corresponding result

%
ind_TH = 10;

%
inds_plt = 1 : datascaling.ind_pairs(ind_TH, 2);

%
datascaling.locID = datascaling.locID(1:ind_TH, :);
%
datascaling.sigmah = datascaling.sigmah(1:ind_TH);
datascaling.ind_pairs = datascaling.ind_pairs(1:ind_TH, :);

%
datascaling.fe = datascaling.fe(inds_plt);
datascaling.Ab = datascaling.Ab(inds_plt);
datascaling.Abonsigmah = datascaling.Abonsigmah(inds_plt);
datascaling.Abonsigmah_ref = datascaling.Abonsigmah_ref(inds_plt);

%
datascaling.r_A_cumulative = datascaling.r_A_cumulative(1:ind_TH);
datascaling.r_B_cumulative = datascaling.r_B_cumulative(1:ind_TH);


%% Do bin-average

%
Abonsigma_grid = logspace(log10(0.2), log10(1), 7);

%
dspacing = diff(Abonsigma_grid);
%
dspacing_vec = [dspacing(1), dspacing];

%
quantiles_range = [0.25, 0.75];

%
fe_mean_sigmah_ref = NaN(size(Abonsigma_grid));
fe_quant_sigmah_ref = NaN(length(Abonsigma_grid), 2);
%
fe_mean_sigmah = NaN(size(Abonsigma_grid));
fe_quant_sigmah = NaN(length(Abonsigma_grid), 2);
%
fe_mean_Ab = NaN(size(Abonsigma_grid));
fe_quant_Ab = NaN(length(Abonsigma_grid), 2);


%
for i = 1:length(Abonsigma_grid)
    
    % --------------------------------------------------
    %
    linbin_aux = (datascaling.Abonsigmah_ref >= (Abonsigma_grid(i) - dspacing_vec(1))) & ...
                 (datascaling.Abonsigmah_ref <= (Abonsigma_grid(i) + dspacing_vec(1)));
             
	%
    fe_mean_sigmah_ref(i) = mean(datascaling.fe(linbin_aux), 'omitnan');
    %
    fe_quant_sigmah_ref(i, :) = quantile(datascaling.fe(linbin_aux), quantiles_range);
    
    
    % --------------------------------------------------
    %
    linbin_aux = (datascaling.Abonsigmah >= (Abonsigma_grid(i) - dspacing_vec(1))) & ...
                 (datascaling.Abonsigmah <= (Abonsigma_grid(i) + dspacing_vec(1)));
             
	%
    fe_mean_sigmah(i) = mean(datascaling.fe(linbin_aux), 'omitnan');
    %
    fe_quant_sigmah(i, :) = quantile(datascaling.fe(linbin_aux), quantiles_range);
    
    % --------------------------------------------------
    %
    linbin_aux = (datascaling.Ab >= (Abonsigma_grid(i) - dspacing_vec(1))) & ...
                 (datascaling.Ab <= (Abonsigma_grid(i) + dspacing_vec(1)));
             
	%
    fe_mean_Ab(i) = mean(datascaling.fe(linbin_aux), 'omitnan');
    %
    fe_quant_Ab(i, :) = quantile(datascaling.fe(linbin_aux), quantiles_range);

end


%% Create file with bin-averaged data and power law-coeffcients

%
feROXSI2022.locationID = datascaling.locID;

%
feROXSI2022.Abonsigma_grid = Abonsigma_grid(:);
%
feROXSI2022.fe_binmean_sigmahvar = fe_mean_sigmah(:);
feROXSI2022.fe_binmean_sigmahref = fe_mean_sigmah_ref(:);

%
feROXSI2022.sigmahvar = datascaling.sigmah;
feROXSI2022.sigmahref = datascaling.sigmah_ref;

%
feROXSI2022.quantiles = quantiles_range;

%
feROXSI2022.fe_quant_sigmahvar = fe_quant_sigmah;
feROXSI2022.fe_quant_sigmahref = fe_quant_sigmah_ref;

%
x_aux = log10(datascaling.Abonsigmah);
y_aux = log10(datascaling.fe);

%
pwlaw_coefs = regress(y_aux, [ones(size(x_aux)), x_aux]);

%
feROXSI2022.coef_1 = 10.^(pwlaw_coefs(1));
feROXSI2022.coef_2 = pwlaw_coefs(2);

%
dir_outputfile = fullfile(paper_directory(), 'data', 'level_5');

%
save(fullfile(dir_outputfile, 'bin_averaged_data.mat'), 'feROXSI2022')


%%
% ----------------------------------------------------------------
% ----------------------------------------------------------------
% ------------------------- MAKE FIGURE --------------------------
% ----------------------------------------------------------------
% ----------------------------------------------------------------

%% Make figure

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.3227    0.1367    0.1886    0.6];
%
haxs = makeSubPlots(0.2, 0.1, 0.1, ...
                    0.03, 0.1, 0.1, 1, 2);
hold(haxs, 'on')

% ---------------------------------------------------------------
% Plot fe vs. Ab/sigma_h

%
mkSZhourly = 8;
mkCLhourly = 0.5.*[1, 1, 1];

    
    %
    plot(haxs(1), datascaling.Abonsigmah(inds_plt), ...
                  datascaling.fe(inds_plt), ...
                  '.', 'MarkerSize', mkSZhourly, 'Color', mkCLhourly)
    %
    plot(haxs(1), Abonsigma_grid, fe_mean_sigmah, '.-k', 'MarkerSize', 28, 'LineWidth', 2)
    

    %
    plot(haxs(2), datascaling.Abonsigmah_ref(inds_plt), ...
                  datascaling.fe(inds_plt), ...
                  '.', 'MarkerSize', mkSZhourly, 'Color', mkCLhourly)
    %
    plot(haxs(2), Abonsigma_grid, fe_mean_sigmah_ref, '.-k', 'MarkerSize', 28, 'LineWidth', 2)
        
%
for i = 1:length(Abonsigma_grid)
    

    % ----------------------------
    %
    dfe_neg_aux = fe_mean_sigmah_ref(i) - fe_quant_sigmah_ref(i, 1);
    dfe_pos_aux = fe_quant_sigmah_ref(i, 2) - fe_mean_sigmah_ref(i);
    %
    herb_aux = errorbar(haxs(2), ...
                        Abonsigma_grid(i), ...
                        fe_mean_sigmah_ref(i), ...
                        dfe_neg_aux, dfe_pos_aux);
        %
        herb_aux.Color = [0, 0, 0];
        herb_aux.LineWidth = 3;
        
    % ----------------------------
    %
    dfe_neg_aux = fe_mean_sigmah(i) - fe_quant_sigmah(i, 1);
    dfe_pos_aux = fe_quant_sigmah(i, 2) - fe_mean_sigmah(i);
    %
    herb_aux = errorbar(haxs(1), ...
                        Abonsigma_grid(i), ...
                        fe_mean_sigmah(i), ...
                        dfe_neg_aux, dfe_pos_aux);
        %
        herb_aux.Color = [0, 0, 0];
        herb_aux.LineWidth = 3;
        

end  


% ------------------------------------------------
% Annotate r2

%
text(haxs(1), 0.11, 0.825, ['$r_*^2 = ' num2str(datascaling.r_A_cumulative(ind_TH).^2, '%.2f') '$'], 'Interpreter', 'Latex', 'FontSize', 20)
text(haxs(2), 0.11, 0.825, ['$r_*^2 = ' num2str(datascaling.r_B_cumulative(ind_TH).^2, '%.2f') '$'], 'Interpreter', 'Latex', 'FontSize', 20)

%
text(haxs(1), 0.104, 15, 'a)', 'Interpreter', 'Latex', 'FontSize', 24)
text(haxs(2), 0.104, 15, 'b)', 'Interpreter', 'Latex', 'FontSize', 24)


% ------------------------------------------------
% Set properties of axes

%
set(haxs, 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
set(haxs, 'XScale', 'log', 'YScale', 'log')
set(haxs, 'XLim', [0.1, 1.5], 'YLim', [7e-1, 20])
%
set(haxs, 'XTick', [1e-2, 1e-1, 0.5, 1e0, 1e1])
set(haxs, 'YTick', [0.1, 1, 5, 10, 100])
   
% -----------------------------------------------
% Add x/y labels

%
hxlb_1 = xlabel(haxs(1), '$\langle A_b \rangle/ \langle \sigma_h \rangle$', 'Interpreter', 'Latex', 'FontSize', 22);
%
hxlb_2 = xlabel(haxs(2), '$\langle A_b \rangle/ \sigma_h^{\mathrm{ref}}$', 'Interpreter', 'Latex', 'FontSize', 22);

%
ylabel(haxs(1), '$f_e$', 'Interpreter', 'Latex', 'FontSize', 22)
ylabel(haxs(2), '$f_e$', 'Interpreter', 'Latex', 'FontSize', 22)

%
hxlb_1.Position(2) = 0.525;
hxlb_2.Position(2) = 0.525;


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', 'figure07.pdf'), 'Resolution', 300)

