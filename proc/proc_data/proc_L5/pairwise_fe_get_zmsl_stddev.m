function pwdissipation = pairwise_fe_get_zmsl_stddev(pwdissipation)
%% pwdissipation = PAIRWISE_FE_GET_ZMSL_STDDEV(pwdissipation)
%
%   inputs
%       - pwdissipation: structure variable created by
%                        pairwise_fe_getdatapair.m
%
%   outputs
%       - pwdissipation: same structure variable, but with the average
%                        depth and standard deviation of depth within a
%                        rectangle bounding the two instrument sites.
%
%
% PAIRWISE_FE_GET_ZMSL_STDDEV.m loads the bathymetry data files and add
% the spatially averaged bottom depth (zmsl) and standard deviation
% (zmsl_stddev_deplaned) within a rectangle bounding instrument sites.
%
% See also:
%   pairwise_fe_getdatapair.m


%% Load bathymetry with roughness

%
dir_bathy = fullfile(paper_directory(), 'data', 'bathymetry');

% Use bathymetry computed with 20 m smoothing (average bottom depth between
% instruments from this file is very similar than when computing from
% the  m gridded bathymetry) 
bathymetry.Asilomar = load(fullfile(dir_bathy, 'bathy_longscale_Asilomar.mat'));
bathymetry.ChinaRock = load(fullfile(dir_bathy, 'bathy_longscale_ChinaRock.mat'));

%
bathymetry.Asilomar = bathymetry.Asilomar.bathymetry;
bathymetry.ChinaRock = bathymetry.ChinaRock.bathymetry;


%% Define extra margin to compute mean
% depth and its variability (roughness)

%
dx = 10;
dy = 10;


%% Loop over pairs of instruments and get bathymetry data

%
for i = 1:pwdissipation.Npairs
    
    %% Get edges of the rectangle bounding instrument pair
    %
    x_ref = mean(pwdissipation.data(i).X, 1, 'omitnan');
    y_ref = mean(pwdissipation.data(i).Y, 1, 'omitnan');
    
    %
    x_edges = [min(x_ref), max(x_ref)] + dx.*[-1, 1];
    y_edges = [min(y_ref), max(y_ref)] + dy.*[-1, 1];
    
    
    %
    lin_xlims_aux = (bathymetry.(pwdissipation.site(i)).x >= x_edges(1)) & ...
                    (bathymetry.(pwdissipation.site(i)).x <= x_edges(2));
    %
    lin_ylims_aux = (bathymetry.(pwdissipation.site(i)).y >= y_edges(1)) & ...
                    (bathymetry.(pwdissipation.site(i)).y <= y_edges(2));

    %
    pwdissipation.data(i).xedges_bathy = x_edges;
    pwdissipation.data(i).yedges_bathy = y_edges;
        
    
    %% Get all zmsl and standard deviation data
    % points within the bounding rectangle
    
    %
    zmsl_aux = bathymetry.(pwdissipation.site(i)).zmsl_mean(lin_ylims_aux, lin_xlims_aux);
    zmsl_aux = zmsl_aux(~isnan(zmsl_aux));

    %
    z_std_aux = bathymetry.(pwdissipation.site(i)).zmsl_stddev_deplaned(lin_ylims_aux, lin_xlims_aux);
    z_std_aux = z_std_aux(~isnan(z_std_aux));
    
    
    %% Average within the bounding rectangle
    
    %
    pwdissipation.data(i).zmsl_all = zmsl_aux;
    pwdissipation.data(i).zmsl_avg = mean(zmsl_aux(:));
    
    %
    pwdissipation.data(i).zmsl_stddev_all = z_std_aux;
    pwdissipation.data(i).zmsl_stddev_avg = mean(z_std_aux(:));
    
    
end