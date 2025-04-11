%% Load data and create info_table.txt for creating Table 1 in latex.

clear
close all


%% Load L5 data

%
dir_data = fullfile(paper_directory(), 'data', 'level_5');

%
load(fullfile(dir_data, 'dataL5_withQC.mat'))


%% Subset table with just sites that passed all QC

%
lpassedallQC = pwfetable.pass_allQC;

%
tablesub = pwfetable(pwfetable.pass_allQC, :);


%% Sort in descending order of the correlation coefficient

%
[~, ind_sort] = sort(tablesub.rcoef, 'descend');

%
tablesub = tablesub(ind_sort, :);


%% Get variables

%
Npairs = size(tablesub, 1);

%
getvars.locationID = [tablesub.locationID_A, tablesub.locationID_B];

%
getvars.hmean = NaN(Npairs, 2);
%
getvars.Hs = NaN(Npairs, 2);
%
getvars.dx = tablesub.delta_X;
getvars.dy = tablesub.delta_Y;
%
getvars.mean_dEcgdx = NaN(Npairs, 1);
getvars.mean_Urms = NaN(Npairs, 2);
%
getvars.mean_Ab = NaN(Npairs, 2);
%
getvars.r2 = tablesub.rcoef.^2;
getvars.bulkfe = tablesub.fe_bulk;
getvars.sigmah = tablesub.sigmah_avg;

%
for i = 1:Npairs
    
    %
    ind_match_A = find(strcmp(pwdissipation.locationID(:, 1), getvars.locationID(i, 1)));
    ind_match_B = find(strcmp(pwdissipation.locationID(:, 2), getvars.locationID(i, 2)));
    
    %
    if ind_match_A~=ind_match_B
        error('!!!!')
    else
        ind_match_aux = ind_match_A;
    end
    
    %
    getvars.hmean(i, :) = pwdissipation.data(ind_match_aux).bottomdepth_mean;
    
    %
    lok_avg = ~isnan(pwdissipation.data(ind_match_aux).fromEcg.dFluxdx);
    
    %
    getvars.Hs(i, :) = mean(pwdissipation.data(ind_match_aux).Hs(lok_avg, :), 1);
    %
    getvars.mean_dEcgdx(i) = mean(pwdissipation.data(ind_match_aux).fromEcg.dFluxdx(lok_avg, :), 1);
    %
    getvars.mean_Urms(i, :) = mean(pwdissipation.data(ind_match_aux).Urms(lok_avg, :), 1);

    %
    getvars.mean_Ab(i, :) = mean(pwdissipation.data(ind_match_aux).Ab(lok_avg, :), 1);
    
end

%
getvars.hmean_avg = (getvars.hmean(:, 1) + getvars.hmean(:, 2))./2;


%% Create *.txt file for saving lines of table
% (If file already exists, it is overwritten)

%
dir_output = fullfile(paper_directory(), 'figures');
%
fid = fopen(fullfile(dir_output, 'info_table.txt'), 'w');



%% Put information together (with average water depth between instruments)

%
fprintf(fid, '%s \n', '\begin{tabular}{cccccccccccc}');
fprintf(fid, '%s \n', '    \hline\hline');
fprintf(fid, '%s \n', ['    $N$ & ID & $\hmean$~[m] & $\Delta x$~[m] & $\Delta y$~[m] & ' ...
                       '$\Hsmean$~[m] & $\overline{U}_{\mathrm{rms}}$~[m s$^{-1}$] & ' ...
                       '$-\overline{\dtxt F / \dtxt x}~[\mathrm{W}\,\mtxt^{-2}]$ & ' ...
                       '$r^2$ & $\bulkf$ & $\overline{A}_b$~[m] & $\mstdh$~[m] \\']);
fprintf(fid, '%s \n', '    \hline');

%
for i = 1:Npairs
    fprintf(fid, ...
            '%s \n', ...
            ['    ' ...
             num2str(i) ' & ' ...
             char(getvars.locationID(i, 1)) '--' char(getvars.locationID(i, 2)) ' & ' ...
             num2str(getvars.hmean_avg(i), '%.1f') ' & ' ...
             num2str(getvars.dx(i), '%.0f') ' & ' num2str(getvars.dy(i), '%.0f')  ' & ' ...
             num2str(num2str(getvars.Hs(i, 1), '%.2f')) '--' num2str(num2str(getvars.Hs(i, 2), '%.2f')) ' & ' ...
             num2str(num2str(getvars.mean_Urms(i, 1), '%.2f')) '--' num2str(num2str(getvars.mean_Urms(i, 2), '%.2f')) ' & ' ...
             num2str(num2str(-getvars.mean_dEcgdx(i, 1), '%.0f')) ' & ' ...
             num2str(num2str(getvars.r2(i), '%.2f')) ' & ' ...
             num2str(num2str(getvars.bulkfe(i), '%.1f')) ' & ' ...
             num2str(num2str(getvars.mean_Ab(i, 1), '%.2f')) '--' num2str(num2str(getvars.mean_Ab(i, 2), '%.2f')) ' & ' ...
             num2str(num2str(getvars.sigmah(i), '%.2f')) ' \\']);
end
fprintf(fid, '%s \n', '    \hline');
fprintf(fid, '%s \n', '\end{tabular}');

%
fclose(fid);
