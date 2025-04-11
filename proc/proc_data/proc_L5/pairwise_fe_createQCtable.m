function fetable = pairwise_fe_createQCtable(pwdissipation, fluxsource)
%% fetable = PAIRWISE_FE_CREATEQCTABLE(pwdissipation, fluxsource)
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
% PAIRWISE_FE_CREATEQCTABLE.m gets data from structure pwdissipation, and
% puts results in a table. Part of PAIRWISE_FE_CREATEQCTABLE.m is the same
% as pairwise_fe_createtable.m, but the Quality Control (QC) table has
% additional variables.
%
%
% See also:
%   run_proc_pairwise_fe_allsites.m
%   pairwise_fe_getdatapair.m
%   pairwise_fe_calculate.m
%   pairwise_fe_createtable.m


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


%% Initialize variables that will be added to the output table

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
h_bathy = NaN(Npairs, 1);

%
sigmah_avg = NaN(Npairs, 1);


%% Loop over pairs

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



%% Create standard table

%
fetable = table(site, locationID_A, locationID_B, ...
                h_A, h_B, h_AB_avg, h_bathy, ...
                X_A, X_B, Y_A, Y_B, X_AB_avg, Y_AB_avg, ...
                delta_X, delta_Y, ratio_dydx, ...
                dFluxdx_mean, fe_bulk, fe_bulk_error, rcoef, ...
                sigmah_avg);

            
%%
% --------------------------------------------------------
% -------- GET VARIABLES SPECIFIC TO THE QC TABLE --------
% --------------------------------------------------------


%% Initialize other variables that will be added to the output table

% ----------------------
%
vec_1_aux = NaN(length(site), 1);
vec_2_aux = false(length(site), 1);

% ----------------------
days_noNaN = vec_1_aux;

% ----------------------
%
fraction_posD = vec_1_aux;
fraction_hdeep = vec_1_aux;
fraction_Hsonhsmall = vec_1_aux;

%
lenough_posD = vec_2_aux;
lenough_hdeep = vec_2_aux;
lenough_Hsonhsmall = vec_2_aux;
%
lmeanposD = vec_2_aux;
lhmeandeepenough = vec_2_aux;
ldxnotsmall = vec_2_aux;
ldxnotlarge = vec_2_aux;
lroughlyonx = vec_2_aux;

% ----------------------
rcoef_dEcgdxvsU3_flagged = vec_1_aux;


%% Fill data into variables initialized above

%
for i = 1:pwdissipation.Npairs
    
    %
    days_noNaN(i) = pwdissipation.data(i).(fluxsource).NptsnoNaN./24;
  
    % ----------------------     
    %
    fraction_posD(i) = pwdissipation.data(i).(fluxsource).fraction_negdFluxdx;
    fraction_hdeep(i) = pwdissipation.data(i).(fluxsource).fraction_hdeep;
    fraction_Hsonhsmall(i) = pwdissipation.data(i).(fluxsource).fraction_Hsonhsmall;
    %
    lenough_posD(i) = pwdissipation.data(i).(fluxsource).lenough_negdFluxdx;
    lenough_hdeep(i) = pwdissipation.data(i).(fluxsource).lenough_hdeep;
    lenough_Hsonhsmall(i) = pwdissipation.data(i).(fluxsource).lenough_Hsonhsmall;

    % ----------------------
           
    %
    lmeanposD(i) = pwdissipation.data(i).(fluxsource).lposDmean;
    lhmeandeepenough(i) = pwdissipation.data(i).lhmeandeepenough;
    %
    ldxnotsmall(i) = pwdissipation.data(i).ldxnotsmall;
    ldxnotlarge(i) = pwdissipation.data(i).ldxnotlarge;
    lroughlyonx(i) = pwdissipation.data(i).lroughlyonx;
                   
    % ----------------------
    %
    rcoef_dEcgdxvsU3_flagged(i) = pwdissipation.data(i).(fluxsource).rcoef;
    
end


%% Create table

%
feQCtable = table(days_noNaN, ...
                  dFluxdx_mean, ...
                  rcoef_dEcgdxvsU3_flagged, ...
                  fraction_posD, fraction_hdeep, fraction_Hsonhsmall, ...
                  lenough_posD, lenough_hdeep, lenough_Hsonhsmall, ...
                  lmeanposD, lhmeandeepenough, ldxnotsmall, ldxnotlarge, lroughlyonx);              


%% Get sites that pall all QC

%
pass_allQC = feQCtable.lenough_posD & ... 
             feQCtable.lenough_hdeep & ...
             feQCtable.lenough_Hsonhsmall & ...
             feQCtable.lmeanposD & ...
             feQCtable.lhmeandeepenough & ...
             feQCtable.ldxnotsmall & ...
             feQCtable.ldxnotlarge & ...
             feQCtable.lroughlyonx;

%
feQCtable.pass_allQC = pass_allQC;


%% Merge 2 tables for the output

%
fetable_fieldremoved = removevars(fetable, 'dFluxdx_mean');

%
fetable = [fetable_fieldremoved, feQCtable];