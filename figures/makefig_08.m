%% Make Figure 8 -- correlations squared as a function of pair number

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

%
ind_plt = 1:size(pwfetable_passedQC, 1);


%% Compute cumulative correlation of log_10(Ab/sigma_h) vs log_10(fe)

%
N_aux = length(ind_passed_QC);

%
sigmah_avgsites = mean(pwfetable_passedQC.sigmah_avg);

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
    
    %% Concatenate data from different pairs in one structure variable
    
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


%% Make figure

%
hfig = figure;
hfig.Units = 'normalized';
hfig.Position = [0.2, 0.2, 0.22, 0.6];
%
haxs = makeSubPlots(0.21, 0.02, 0.1, ...
                    0.1, 0.1, 0.05, 1, 2);
hold(haxs, 'on')
 
    %
    plot(haxs(1), ind_plt, pwfetable_passedQC.rcoef.^2, '.-k', 'MarkerSize', 28, 'LineWidth', 2)
    
    %
    indtrimplt = 2:length(ind_plt);
    %
    plot(haxs(2), ind_plt(indtrimplt), datascaling.r_A_cumulative(indtrimplt).^2, '.-', 'MarkerSize', 28, 'LineWidth', 2)
    plot(haxs(2), ind_plt(indtrimplt), datascaling.r_B_cumulative(indtrimplt).^2, '.-', 'MarkerSize', 28, 'LineWidth', 2)
    
    
%
hleg = legend(haxs(2), '$\langle \sigma_h \rangle$', ...
                       '$\sigma_h^{\mathrm{ref}}$', ...
                       'Location', 'SouthEast', 'Interpreter', 'Latex', 'FontSize', 22);
hleg.Position = [0.7671    0.1047    0.2052    0.1054];
    

% --------------------------------------------
% Set axes

%
set(haxs, 'FontSize', 14, 'Box', 'on', 'XGrid', 'on', 'YGrid', 'on')
set(haxs, 'XLim', [0.5, (ind_plt(end)+0.5)])
set(haxs, 'XTick', 1:length(ind_plt))
%
set(haxs(1), 'XTickLabel', [])
%
set(haxs(1), 'YLim', [0, 1])
set(haxs(2), 'YLim', [0, 0.5])

%
hxaxs = get(haxs(2), 'XAxis');
hxaxs.FontSize = 14;
hxaxs.TickLabelRotation = 45;

% --------------------------------------------
% Add x/y labels

%
text(haxs(1), 0.7, 0.07, 'a)', 'Interpreter', 'Latex', 'FontSize', 24)
text(haxs(2), 0.7, 0.46, 'b)', 'Interpreter', 'Latex', 'FontSize', 24)

%
xlabel(haxs(2), 'Pair number $N$', 'Interpreter', 'Latex', 'FontSize', 18)
%
ylabel(haxs(1), '$r^2 (N)$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel(haxs(2), '$r_*^2 (1$:$N)$', 'Interpreter', 'Latex', 'FontSize', 18)

%
linkaxes(haxs, 'x')


%% Save figure

%
exportgraphics(hfig, fullfile(paper_directory(), 'figures', ...
                              'figure08.pdf'), 'Resolution', 300)


