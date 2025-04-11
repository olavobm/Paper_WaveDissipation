%% High-legel script that calls all scripts to make all figures
% (and generates text for table in latex).

clear
close all


%% Add directory/subdirectories to path

%
addpath(genpath(paper_directory()))


%% Figure 01 -- maps

%
run(fullfile(paper_directory(), 'figures', 'makefig_01.m'))


%% Figure 02 -- standard deviation of bathymetry

%
run(fullfile(paper_directory(), 'figures', 'makefig_02.m'))


%% Figure 03 -- timeseries with overview of wave conditions

%
run(fullfile(paper_directory(), 'figures', 'makefig_03.m'))


%% Figure 04 -- cross-shore profiles of statistics

%
run(fullfile(paper_directory(), 'figures', 'makefig_04.m'))


%% Figure 05 -- timeseries example with bulk fe (and map)

%
run(fullfile(paper_directory(), 'figures', 'makefig_05.m'))


%% Figure 06 -- sigma_h corresponding to Fig. 5

%
run(fullfile(paper_directory(), 'figures', 'makefig_06.m'))


%% Table 1 -- run script that gets data for table and create text that
% can be copied and pasted into Latex

run(fullfile(paper_directory(), 'figures', 'make_table_results.m'))


%% Figure 07 -- hourly fe vs Ab/sigmah

%
run(fullfile(paper_directory(), 'figures', 'makefig_07.m'))


%% Figure 08 -- sensitivity of r_*^2

run(fullfile(paper_directory(), 'figures', 'makefig_08.m'))


%% Figure 09 -- dF/dx vs. dFx/dx

run(fullfile(paper_directory(), 'figures', 'makefig_09.m'))


%% Figure 10 -- plot figure with theoretical parameterizations and data

run(fullfile(paper_directory(), 'figures', 'makefig_10.m'))

