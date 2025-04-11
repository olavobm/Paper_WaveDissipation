function fetable = pairwise_fe_createtable(pwdissipation, fluxsource)
%% fetable = PAIRWISE_FE_CREATETABLE(pwdissipation, fluxsource)
%
%   inputs
%       - pwdissipation: data structure by script
%                        run_proc_pairwise_fe_allsites.m
%       - fluxsource (optional): string specifying which flux variable will
%                                be used to compute flux convergence
%                                ("fromEcg" or "fromFx").
%
%   outputs
%       - fetable: a table with fe results.
%
%
% PAIRWISE_FE_CREATETABLE.m gets data from structure pwdissipation, and
% puts results in a table.
%
%
% See also:
%   run_proc_pairwise_fe_allsites.m
%   pairwise_fe_getdatapair.m
%   pairwise_fe_calculate.m


%% If fluxsource is not given, use fromEcg

if ~exist('fluxsource', 'var')
    fluxsource = "fromEcg";
end


%% Get location IDs

%
site = pwdissipation.site;

%
locationID_A = pwdissipation.locationID(:, 1);
locationID_B = pwdissipation.locationID(:, 2);


%% Initialize all variables that will be added to the output table

%
Npairs = pwdissipation.Npairs;

%
h_A = NaN(Npairs, 1);
h_B = NaN(Npairs, 1);

%
X_A = NaN(Npairs, 1);
X_B = NaN(Npairs, 1);
%
Y_A = NaN(Npairs, 1);
Y_B = NaN(Npairs, 1);

%
Hs_A = NaN(Npairs, 1);     Hs_B = NaN(Npairs, 1);
Urms_A = NaN(Npairs, 1);   Urms_B = NaN(Npairs, 1);
Ab_A = NaN(Npairs, 1);     Ab_B = NaN(Npairs, 1);

%
Flux_A = NaN(Npairs, 1);    Flux_B = NaN(Npairs, 1);
%
dFluxdx_mean = NaN(Npairs, 1);

%
fe_bulk = NaN(Npairs, 1);

%
fe_bulk_error = NaN(Npairs, 1);

%
rcoef = NaN(Npairs, 1);


%
xedges_bathy = NaN(Npairs, 2);
yedges_bathy = NaN(Npairs, 2);

%
h_bathy = NaN(Npairs, 1);

%
sigmah_avg = NaN(Npairs, 1);


%% Loop over instruments pairs

%
for i = 1:Npairs
    
    %
    lok_common = ~isnan(pwdissipation.data(i).(fluxsource).Flux(:, 1)) & ...
                 ~isnan(pwdissipation.data(i).(fluxsource).Flux(:, 2));
    
    %
    h_A(i) = mean(pwdissipation.data(i).bottomdepth(lok_common, 1));
    h_B(i) = mean(pwdissipation.data(i).bottomdepth(lok_common, 2));
    
    %
    X_A(i) = mean(pwdissipation.data(i).X(lok_common, 1));
    X_B(i) = mean(pwdissipation.data(i).X(lok_common, 2));
    %
    Y_A(i) = mean(pwdissipation.data(i).Y(lok_common, 1));
    Y_B(i) = mean(pwdissipation.data(i).Y(lok_common, 2));

    %
    Hs_A(i) = mean(pwdissipation.data(i).Hs(lok_common, 1));
    Hs_B(i) = mean(pwdissipation.data(i).Hs(lok_common, 2));
    %
    Urms_A(i) = mean(pwdissipation.data(i).Urms(lok_common, 1));
    Urms_B(i) = mean(pwdissipation.data(i).Urms(lok_common, 2));
    %
    Ab_A(i) = mean(pwdissipation.data(i).Ab(lok_common, 1));
    Ab_B(i) = mean(pwdissipation.data(i).Ab(lok_common, 2));
    %
    Flux_A(i) = mean(pwdissipation.data(i).(fluxsource).Flux(lok_common, 1));
    Flux_B(i) = mean(pwdissipation.data(i).(fluxsource).Flux(lok_common, 2));
    %
    dFluxdx_mean(i) = mean(pwdissipation.data(i).(fluxsource).dFluxdx(lok_common));
    
    
    %% Get bulk fe and associated quantities

    %
    fe_bulk(i) = pwdissipation.data(i).(fluxsource).posD.fctefit;

    %
    fe_bulk_error(i) = pwdissipation.data(i).(fluxsource).posD.fcteerror(2) - ...
                       pwdissipation.data(i).(fluxsource).posD.fctefit;

    %
    rcoef(i) = pwdissipation.data(i).(fluxsource).posD.rcoef;
    
    
    %% Get variables related to bathymetry
    
    %
    xedges_bathy(i, :) = pwdissipation.data(i).xedges_bathy;
    yedges_bathy(i, :) = pwdissipation.data(i).yedges_bathy;
    
    % 
    h_bathy(i) = -pwdissipation.data(i).zmsl_avg;

    %
    sigmah_avg(i) = pwdissipation.data(i).zmsl_stddev_avg;
    
end


%% Calculate quantities based on the variables retrieved in the loop

%
h_AB_avg = (h_A + h_B)./2;

%
X_AB_avg = (X_A + X_B)./2;
Y_AB_avg = (Y_A + Y_B)./2;

%
delta_X = X_B - X_A;
delta_Y = Y_B - Y_A;

%
ratio_dydx = delta_Y ./ delta_X;


%% Create table

fetable = table(site, locationID_A, locationID_B, ...
                h_A, h_B, h_AB_avg, h_bathy, ...
                X_A, X_B, Y_A, Y_B, X_AB_avg, Y_AB_avg, ...
                delta_X, delta_Y, ratio_dydx, ...
                dFluxdx_mean, fe_bulk, fe_bulk_error, rcoef, ...
                xedges_bathy, yedges_bathy, sigmah_avg);