%% Compute standard deviation from point clound bathymetry
% (and also compute a mean depth using this long averaging scale)

clear
close all



%%  Set directories

%
dir_data = fullfile(paper_directory(), 'data', 'bathymetry');

%
dir_output = dir_data;


%% Parameters for China Rock and Asilomar

%
list_sites = ["ChinaRock", "Asilomar"];


% ---------------------------------------------
%
gridinfo.ChinaRock.filename_pointcloud = ["zmsl_combined_ChinaRock.mat"];
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
gridinfo.Asilomar.filename_pointcloud = ["zmsl_combined_Asilomar.mat"];
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
disp(['Start processing at: ' datestr(now())])


%% Now compute mean and deplaned standard deviation within the grid cells

%
for i1 = 1:length(list_sites)
    
    %% Load point cloud
    
    elevroxsi = load(fullfile(paper_directory(), 'data', 'bathymetry', ...
                              gridinfo.(list_sites(i1)).filename_pointcloud));
	elevroxsi = elevroxsi.elevroxsi;


    %% Start structure array
    
    %
    bathymetry.site = list_sites(i1);
    
    %
    bathymetry.xlims = gridinfo.(list_sites(i1)).xedges;
    bathymetry.ylims = gridinfo.(list_sites(i1)).yedges;
    
    %
    bathymetry.dx = gridinfo.(list_sites(i1)).dx;
    %
    bathymetry.dx_window = gridinfo.(list_sites(i1)).dx_window;
    %
    bathymetry.Npts_TH = gridinfo.(list_sites(i1)).Npts_TH;
    
    %
    bathymetry.x = gridinfo.(list_sites(i1)).x;
    bathymetry.y = gridinfo.(list_sites(i1)).y;
    
    %
    prealloc_aux = NaN(length(bathymetry.y), length(bathymetry.x));
    
    %
    bathymetry.Npts = prealloc_aux;
    
    %
    bathymetry.zmsl_mean = prealloc_aux;
    bathymetry.zmsl_stddev = prealloc_aux;
    %
    bathymetry.zmsl = prealloc_aux;
    bathymetry.zmsl_stddev_deplaned = prealloc_aux;
   
    
    %% Remove from point cloud points that are far away from the grid
    
    %
    lkeep_aux = (elevroxsi.x >= (bathymetry.xlims(1) - bathymetry.dx_window)) & ...
                (elevroxsi.x <= (bathymetry.xlims(2) + bathymetry.dx_window)) & ...
                (elevroxsi.y >= (bathymetry.ylims(1) - bathymetry.dx_window)) & ...
                (elevroxsi.y <= (bathymetry.ylims(2) + bathymetry.dx_window));
    
    %
    list_fields_aux = ["latitude", "longitude", "easting", "northing", "x", "y", "z_msl"];
    %
    for i_aux = 1:length(list_fields_aux)
        %
        elevroxsi.(list_fields_aux(i_aux)) = elevroxsi.(list_fields_aux(i_aux))(lkeep_aux);
    end
            
            
    %% Loop over grid points, compute mean, (deplaned) standard deviation,
    % and get statistical quantities
    
    %
    for i2 = 1:length(gridinfo.(list_sites(i1)).y)
        
        %
        lin_y_aux = (elevroxsi.y >= (gridinfo.(list_sites(i1)).y(i2) - (gridinfo.(list_sites(i1)).dx_window/2))) & ...
                    (elevroxsi.y <= (gridinfo.(list_sites(i1)).y(i2) + (gridinfo.(list_sites(i1)).dx_window/2)));
        
        %
        x_sub_aux = elevroxsi.x(lin_y_aux);
        y_sub_aux = elevroxsi.y(lin_y_aux);
        z_sub_aux = elevroxsi.z_msl(lin_y_aux);
                
        %
        for i3 = 1:length(gridinfo.(list_sites(i1)).x)
            
            %
            lin_x_aux = (x_sub_aux >= (gridinfo.(list_sites(i1)).x(i3) - (gridinfo.(list_sites(i1)).dx_window/2))) & ...
                        (x_sub_aux <= (gridinfo.(list_sites(i1)).x(i3) + (gridinfo.(list_sites(i1)).dx_window/2)));
        
            %
            x_incell_aux = x_sub_aux(lin_x_aux);
            y_incell_aux = y_sub_aux(lin_x_aux);
            z_incell_aux = z_sub_aux(lin_x_aux);
            
            %
            Npts_aux = length(z_incell_aux);
            
            
            %%

            %
            bathymetry.Npts(i2, i3) = Npts_aux;
            
            %
            if Npts_aux < gridinfo.Asilomar.Npts_TH
                continue 
            end
           
            
            %% Compute plane trend
            
            %
            x_anom_aux = x_incell_aux - gridinfo.(list_sites(i1)).x(i3);
            %
            y_anom_aux = y_incell_aux - gridinfo.(list_sites(i1)).y(i2);
            
            %
            planecoefs_aux = regress(z_incell_aux, [ones(Npts_aux, 1), x_anom_aux, y_anom_aux]);
            
            %
            z_plane_aux = planecoefs_aux(1) + planecoefs_aux(2).*x_anom_aux + ...
                                              planecoefs_aux(3).*y_anom_aux;
            
            %
            z_anom_aux = z_incell_aux - z_plane_aux;
            
            
            %% Compute mean and deplaned standard deviation
           
            %
            bathymetry.zmsl_mean(i2, i3) = mean(z_incell_aux);
            bathymetry.zmsl_stddev(i2, i3) = std(z_incell_aux);
            
            %
            bathymetry.zmsl(i2, i3) = planecoefs_aux(1);      
            
            %
            bathymetry.zmsl_stddev_deplaned(i2, i3) = sqrt(mean(z_anom_aux.^2));
            
        end
    end
 
    
    %% Save file
    
    save(fullfile(dir_output, ['bathy_longscale_' char(list_sites(i1)) '.mat']), 'bathymetry', '-v7.3')
    
    
    %%
    
    clear bathy
    
    %% Print progress message
    
    %
    disp(['done with ' char(list_sites(i1)) ' at: ' char(datetime("now"))])
    
end








