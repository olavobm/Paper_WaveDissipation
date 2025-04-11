%% Script that calls all scripts that do the data processing for the paper.

clear
close all


%% Add directory/subdirectories to path

%
addpath(genpath(paper_directory()))


%%
% ----------------------------------------------
% ---------- L1-to-L2 DATA PROCESSING ----------
% ----------------------------------------------

%% SoloD

%
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L2', 'run_proc_SoloD_L2.m'))


%% Spotter (includes SmartMoorings)

%
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L2', 'run_proc_Spotter_L2.m'))


%% ADCPs (Aquadopp and Signature) -- L1-to-L2 processing

%
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L2', 'run_proc_ADCP_L2.m'))


%%
% ----------------------------------------------
% ---------- L2-to-L3 DATA PROCESSING ----------
% ----------------------------------------------

%% Export L2 to L3 data file(s)

%
disp(' '), disp(' ')
%
disp('# Exporting data from L2 to L3')

%
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L3', 'export_L2toL3.m'))
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L3', 'merge_L3_files.m'))

%
disp('--- Done putting all data in a single structure ---')


%% Effective depth correction

%
disp('# Computing and applying effective depth correction')

%
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L3', 'convert_L3data_hp_to_heff.m'))

%
disp('--- Done with effective depth correction ---')

%%
% ----------------------------------------------
% ---------- L3-to-L4 DATA PROCESSING ----------
% ----------------------------------------------


%% L3 to L4


%
disp('# Computing sea-swell bulk statistics')

%
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L4', 'compute_Suu_SAb.m'))
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L4', 'compute_bulkstatistics.m'))

% Do the same for files with the "observed" depth instead of heff
% (this is exactly the same as the 2 scripts above, the only
% difference is that it points to a different L3 file and saves a L4
% file with another name)
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L4', 'proc_hpdata_L3toL4.m'))

%
disp('--- Done with sea-swell bulk statistics ---')


%%
%
disp(['--- DONE WITH LOW-LEVEL PROCESSING AT: ' char(string(datetime('now'))) ' ---'])
               


%%
% ----------------------------------------------------------------------
% -------------------- HIGHER-LEVEL DATA PROCESSING --------------------
% ----------------------------------------------------------------------


%% Compute pair-wise fe for all locations

%
run(fullfile(paper_directory(), 'proc', 'proc_data', 'proc_L5', ...
                                'run_proc_pairwise_fe_allsites.m'))


%% Compute and add bathymetric roughness to pairwise fe

%
pwdissipation = pairwise_fe_get_zmsl_stddev(pwdissipation);


%% Put all results in a table

%
pwfetable = pairwise_fe_createtable(pwdissipation);


%% Re-do fe analysis applying QC criteria and create a QC table
% with results. Also create a table analogous to the one created above

%
pwdissipation_flagged = pairwise_fe_QCflag(pwdissipation);

%
pwdissipation_flagged = pairwise_fe_calculate(pwdissipation_flagged);

%
pwfeQCtable_flagged = pairwise_fe_createQCtable(pwdissipation_flagged);


%% Save L5 files

%
save(fullfile(paper_directory(), 'data', 'level_5', 'dataL5_noQC.mat'), ...
     'pwdissipation', 'pwfetable')
 
%
pwdissipation = pwdissipation_flagged;
pwfetable = pwfeQCtable_flagged;
 
%
save(fullfile(paper_directory(), 'data', 'level_5', 'dataL5_withQC.mat'), ...
     'pwdissipation', 'pwfetable')


%% Print message that data processing is done

%
disp(['--- DONE WITH HIGH-LEVEL PROCESSING AT: ' char(string(datetime('now'))) ' ---'])
%
disp(' ')
%
disp(['####### DONE WITH ALL DATA PROCESSING AT: ' char(string(datetime('now'))) ' #######'])
               



