%% Run this script to reproduce all analyses and figures.
%
% High-legel script that calls other high-level scripts
% for doing all analses and making figures.
%
% Olavo Badaro Marques

clear
close all

%% Add directory/subdirectories to path

%
addpath(genpath(paper_directory()))

%% 1) Process bathymetry data

%
run(fullfile(paper_directory(), 'proc', 'proc_bathymetry', ...
                                'compute_2mscale_from_pointcloud.m'))
%
run(fullfile(paper_directory(), 'proc', 'proc_bathymetry', ...
                                'compute_20mscale_from_pointcloud.m'))
%
run(fullfile(paper_directory(), 'proc', 'proc_bathymetry', ...
                                'compute_heff_maps.m'))

%% 2) Process data from instruments deployed
% in the large-scale array of ROXSI2022

%
run(fullfile(paper_directory(), 'proc', 'proc_data', ...
                                'run_all_dataproc.m'))

%% 3) Make figures

%
run(fullfile(paper_directory(), 'figures', 'run_all_makefigures.m'))

%%
char(datetime('now'))