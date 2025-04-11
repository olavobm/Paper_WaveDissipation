%% Calculate a smoothed bathymetry (effective depth) and
% a map of depth correction for the instruments.

%
clear
close all


%%  Set directories

%
dir_data = fullfile(paper_directory(), 'data', 'bathymetry');

%
dir_output = dir_data;


%% Parameters for computing effective depth.
%
% Parameters are based on results from Marques et al. (2024),
% published at JTECH (doi.org/10.1175/JTECH-D-23-0118.1).
%
% In this paper, effective depth was computed for observations around
% h = 10 m bottom depth, based on measurements of significant wave height
% (Hs). The upper-bound frequency for computing was 0.2 Hz. At this
% frequency and h = 10 m, the wavelength computed from linear-wave theory
% is approximately lambda = 37 m. The optimal radius for averaging the
% bathymetry was found in the paper to be r = 13 m.
%
% For obtaining an averaging radius r for bottom depths different tha
% 10 m, lambda is computed for other bottom depths, and the ratio r/lambda
% is kept constant (r_to_lambda_cte). 

%
freqref = 0.2;
r_to_lambda_cte = 13/37;


%% Parameters for China Rock and Asilomar

%
list_sites = ["ChinaRock", "Asilomar"];

% ---------------------------------------------
%
gridinfo.ChinaRock.filename_pointcloud = "zmsl_combined_ChinaRock.mat";
gridinfo.ChinaRock.filename_longscale = "bathy_longscale_ChinaRock.mat";
gridinfo.ChinaRock.filename_highres = "bathy_2m_ChinaRock.mat";
%
gridinfo.ChinaRock.dx = 2;
%
gridinfo.ChinaRock.dx_window = 20;
%
gridinfo.ChinaRock.Npts_TH = 20;

%
gridinfo.ChinaRock.xedges = [-1000, 20];
gridinfo.ChinaRock.yedges = [-600, 600];

% ---------------------------------------------
%
gridinfo.Asilomar.filename_pointcloud = "zmsl_combined_Asilomar.mat";
gridinfo.Asilomar.filename_longscale = "bathy_longscale_Asilomar.mat";
gridinfo.Asilomar.filename_highres = "bathy_2m_Asilomar.mat";
%
gridinfo.Asilomar.dx = 2;
%
gridinfo.Asilomar.dx_window = 20;
%
gridinfo.Asilomar.Npts_TH = 20;

%
gridinfo.Asilomar.xedges = [-1000, 20];
gridinfo.Asilomar.yedges = [-150, 150];

% ---------------------------------------------
%
for i = 1:length(list_sites)
    
    %
    gridinfo.(list_sites(i)).x = gridinfo.(list_sites(i)).xedges(1) : ...
                                 gridinfo.(list_sites(i)).dx : ...
                                 gridinfo.(list_sites(i)).xedges(2);
    %
    gridinfo.(list_sites(i)).y = gridinfo.(list_sites(i)).yedges(1) : ...
                                 gridinfo.(list_sites(i)).dx : ...
                                 gridinfo.(list_sites(i)).yedges(2);
                             
	%
    gridinfo.(list_sites(i)).x = gridinfo.(list_sites(i)).x(:);
    gridinfo.(list_sites(i)).y = gridinfo.(list_sites(i)).y(:);
                             
end


%
disp(['Start processing at at: ' datestr(now())])


%%
% --------------------------------------------------------------
% --------------------------------------------------------------
% --------------------------------------------------------------

%%

%
for i1 = 1:length(list_sites)
    
    %% Load bathymetry files

    % ---------------------------------
    %
    bathypointcloud = load(fullfile(paper_directory(), ...
                                    'data', 'bathymetry', ...
                                    gridinfo.(list_sites(i1)).filename_pointcloud));
    bathypointcloud = bathypointcloud.elevroxsi;

    % ---------------------------------
    %
    bathyhighres = load(fullfile(paper_directory(), ...
                                 'data', 'bathymetry', ...
                                 gridinfo.(list_sites(i1)).filename_highres));
    bathyhighres = bathyhighres.bathymetry;
    
    
    % ---------------------------------
    %
    bathysgriddedlong = load(fullfile(paper_directory(), ...
                                      'data', 'bathymetry', ...
                                      gridinfo.(list_sites(i1)).filename_longscale));
    bathysgriddedlong = bathysgriddedlong.bathymetry;
    
    
    %% Calculate a grid of the averaging radius
    
    % Don't do the calculation at very shallow depths,
    % where the wavenumber calculation will fail (and where
    % the result is not applicable to the available measurements)
    z_msl_TH = -0.5;
    %
    lbadTH = (bathysgriddedlong.zmsl > z_msl_TH);
    
    %
    k_all = NaN(size(bathysgriddedlong.zmsl));
    
    % Compute wavelength
    for i2 = 1:length(bathysgriddedlong.y)        
        %
        for i3 = 1:length(bathysgriddedlong.x)

            %
            if lbadTH(i2, i3)
                continue
            end
            
            %
            k_all(i2, i3) = wave_freqtok(freqref, -bathysgriddedlong.zmsl(i2, i3));
            
        end 
    end
    
    %
    lambda_all = 2*pi./k_all;
    
    %
    r_map = r_to_lambda_cte .* lambda_all;
    
    
    %% Assign variables to new data structure

    %
    bathymetry.site = list_sites(i1);
    
    %
    bathymetry.dx = bathyhighres.dx;
    %
    bathymetry.x = bathyhighres.x;
    bathymetry.y = bathyhighres.y;
     
    %
    Nx_aux = length(bathymetry.x);
    Ny_aux = length(bathymetry.y);
    
    %
    bathymetry.Npts = NaN(Ny_aux, Nx_aux);
    %
    bathymetry.r = r_map;
    
    %
    bathymetry.Npts_inr = NaN(Ny_aux, Nx_aux);
    %
    bathymetry.zmsl = NaN(Ny_aux, Nx_aux);
    

    %% First get rid of data in pointcloud that are far away from the grid

    %
    rTH = 1.1.*max(r_map(:));
    
    %
    linarea = (bathypointcloud.x >= (bathymetry.x(1) - rTH)) & ...
              (bathypointcloud.x <= (bathymetry.x(end) + rTH)) & ...
              (bathypointcloud.y >= (bathymetry.y(1) - rTH)) & ...
              (bathypointcloud.y <= (bathymetry.y(end) + rTH));
    
	%
    list_aux = ["latitude", "longitude", "easting", "northing", ...
                "x", "y", "z_msl"];
	%
    for i2 = 1:length(list_aux)
        %
        bathypointcloud.(list_aux(i2)) = bathypointcloud.(list_aux(i2))(linarea);
    end
          

    %% Smooth bathymetry with varying radius

    %
    ypts_progress_msg = round(linspace(1, Ny_aux, 10));

    %
    for i2 = 1:Ny_aux

        %
        max_r_i2 = max(bathymetry.r(i2, :));

        %
        liny_aux = (bathypointcloud.y >= (bathymetry.y(i2) - max_r_i2)) & ...
                   (bathypointcloud.y <= (bathymetry.y(i2) + max_r_i2));

        %
        x_aux = bathypointcloud.x(liny_aux);
        y_aux = bathypointcloud.y(liny_aux);
        z_aux = bathypointcloud.z_msl(liny_aux);

        %
        for i3 = 1:Nx_aux

            %
            dist_aux = sqrt((x_aux - bathymetry.x(i3)).^2 + ...
                            (y_aux - bathymetry.y(i2)).^2);

            %
            lindist_aux = (dist_aux <= r_map(i2, i3));

            %
            xinr_aux = x_aux(lindist_aux);
            yinr_aux = y_aux(lindist_aux);
            zinr_aux = z_aux(lindist_aux);

            %
            bathymetry.Npts_inr(i2, i3) = length(find(lindist_aux));

            %
            if bathymetry.Npts_inr(i2, i3) >= 30

                % Use median (as in Marques et al. (2024) JTECH paper)
                bathymetry.zmsl(i2, i3) = median(zinr_aux);
            end

        end

        %
        if any(ypts_progress_msg == i2)
            disp(['Done with gridding data along y grid line ' num2str(i2) ' out ' ...
                  'of ' num2str(Ny_aux) ', at ' datestr(now()) '.'])
        end
        
    end
    
    
    %% Now get the effective depth correction: a factor that needs
    % to be added to the observed (positive) bottom depth measured
    % by pressure sensors.

    %
    bathymetry.dh_bathy = (-bathymetry.zmsl) - (-bathyhighres.zmsl);

    % Remove some very bad points
    bathymetry.dh_bathy(abs(bathymetry.dh_bathy) > 6) = NaN;
    
    
    %% Save file
    
    %
    file_output = fullfile(dir_output, ['bathy_dhbathy_' char(list_sites(i1)) '.mat']);
    
    %
    save(file_output, 'bathymetry')
    
    %
    disp(['--- Done with effective bathymetry at ' char(list_sites(i1)) ' ---'])

    %
    clear bathymetry
 
end

% Print final progress message
disp(['Done with smoothed bathymetry for effective depth at: ' datestr(now())])




