function pwdissipation = pairwise_fe_QCflag(pwdissipation, fluxsource)
%% pwdissipation = PAIRWISE_FE_QCFLAG(pwdissipation, fluxsource)
%
%   inputs
%       - pwdissipation: data structure by script
%                        run_proc_pairwise_fe_allsites.m
%       - fluxsource (optional): string specifying which flux variable will
%                                be used to compute flux convergence
%                                ("fromEcg" or "fromFx").
%
%   outputs
%       - pwdissipation: similar as input, but with QC-flagged results.
%
%
% PAIRWISE_FE_QCFLAG.m flags results in the data structure pwdissipation
% according to several Quality Control (QC) criteria defined below.


%% If fluxsource is not given, use fromEcg

if ~exist('fluxsource', 'var')
    fluxsource = "fromEcg";
end


%% Get data points that are not NaN (in data from both locations
% in each pair) before applying any QC criteria

%
for i = 1:pwdissipation.Npairs

    %
    pwdissipation.data(i).(fluxsource).lnoNaNbeforeQC = ...
                              ~isnan(pwdissipation.data(i).Urms_avg) & ...
                              ~isnan(pwdissipation.data(i).(fluxsource).dFluxdx);
    %
    pwdissipation.data(i).(fluxsource).NptsnoNaN = ...
                 length(find(pwdissipation.data(i).(fluxsource).lnoNaNbeforeQC));
    
end


%% Define QC criteria

%
QCcriteria.dxmin = 20;
QCcriteria.dxmax = 120;
%
QCcriteria.dyondxmax = tand(30);
%
QCcriteria.hmin = 2;
%
QCcriteria.Hsonhmax = 0.25;
%
QCcriteria.r_U3dFxdx = 0.7;

%
QCcriteria.fraction_posD = 0.8;
QCcriteria.fraction_hmin = 0.8;
QCcriteria.fraction_Hsonh = 0.8;

% Add to data structure
pwdissipation.QCcriteria = QCcriteria;


%% Calculate parameters for QC timeseries criteria

%
for i = 1:pwdissipation.Npairs
   
    %% Compute ratio between Hs and bottom depth h
    
    pwdissipation.data(i).Hsonh = pwdissipation.data(i).Hs ./ ...
                                  pwdissipation.data(i).bottomdepth;
    

    %% Get logical variables applying QC criteria to timeseries
    % (i.e., some fraction of the timeseries is selected as passing QC)
    
    %
    lposD = pwdissipation.data(i).(fluxsource).lposD;
    
    %
    ldeepenough = (pwdissipation.data(i).bottomdepth(:, 1) > QCcriteria.hmin) & ...
                  (pwdissipation.data(i).bottomdepth(:, 2) > QCcriteria.hmin);
    %
    lnotbreaking = (pwdissipation.data(i).Hsonh(:, 1) < QCcriteria.Hsonhmax) & ...
                   (pwdissipation.data(i).Hsonh(:, 2) < QCcriteria.Hsonhmax);
    
	%
    pwdissipation.data(i).(fluxsource).lnegdFluxdx = lposD;
    pwdissipation.data(i).ldeepenough = ldeepenough;
    pwdissipation.data(i).lnotbreaking = lnotbreaking;


    %% Calculate the fractions of the logical variables determined above

    %
    pwdissipation.data(i).(fluxsource).fraction_negdFluxdx = ...
                    length(find(pwdissipation.data(i).(fluxsource).lnoNaNbeforeQC & ...
                                pwdissipation.data(i).(fluxsource).lnegdFluxdx)) ./ ...
                    pwdissipation.data(i).(fluxsource).NptsnoNaN;
	%
    pwdissipation.data(i).(fluxsource).fraction_hdeep = ...
                        length(find(pwdissipation.data(i).(fluxsource).lnoNaNbeforeQC & ...
                                    pwdissipation.data(i).ldeepenough)) ./ ...
                        pwdissipation.data(i).(fluxsource).NptsnoNaN;
	%
    pwdissipation.data(i).(fluxsource).fraction_Hsonhsmall = ...
                    length(find(pwdissipation.data(i).(fluxsource).lnoNaNbeforeQC & ...
                                pwdissipation.data(i).lnotbreaking)) ./ ...
                    pwdissipation.data(i).(fluxsource).NptsnoNaN;

    %
    fraction_negdFluxdx = pwdissipation.data(i).(fluxsource).fraction_negdFluxdx;
    fraction_hdeep = pwdissipation.data(i).(fluxsource).fraction_hdeep;
    fraction_Hsonhsmall = pwdissipation.data(i).(fluxsource).fraction_Hsonhsmall;

	%
    pwdissipation.data(i).(fluxsource).fraction_negdFluxdx = fraction_negdFluxdx;
    pwdissipation.data(i).(fluxsource).fraction_hdeep = fraction_hdeep;
    pwdissipation.data(i).(fluxsource).fraction_Hsonhsmall = fraction_Hsonhsmall;
	%
    pwdissipation.data(i).(fluxsource).lenough_negdFluxdx = (fraction_negdFluxdx >= QCcriteria.fraction_posD);
    pwdissipation.data(i).(fluxsource).lenough_hdeep = (fraction_hdeep >= QCcriteria.fraction_hmin);
    pwdissipation.data(i).(fluxsource).lenough_Hsonhsmall = (fraction_Hsonhsmall >= QCcriteria.fraction_Hsonh);


	%% Combine all criteria above and determine which data points
    % pass all QC criteria

    %
    pwdissipation.data(i).(fluxsource).lalltmsrsQC = lposD & ldeepenough & lnotbreaking;
    
    
end


%% Calculate parameters for QC single-parameter criteria

% Loop over pairs of instruments
for i = 1:pwdissipation.Npairs
   
    %
    pwdissipation.data(i).(fluxsource).dFluxdx_mean = ...
                mean(pwdissipation.data(i).(fluxsource).dFluxdx, 'omitnan');
    
    %
    pwdissipation.data(i).(fluxsource).lposDmean = (pwdissipation.data(i).(fluxsource).dFluxdx_mean < 0);
    
    %
    x_mean = mean(pwdissipation.data(i).X, 1, 'omitnan');
    y_mean = mean(pwdissipation.data(i).Y, 1, 'omitnan');
    
    %
    delta_x_mean = diff(x_mean);
    delta_y_mean = diff(y_mean);
    %
    dyondx = delta_y_mean ./ delta_x_mean;
    
    % Get which pairs are not spaced too closely, and not too far apart
    pwdissipation.data(i).ldxnotsmall = (delta_x_mean > QCcriteria.dxmin);
    pwdissipation.data(i).ldxnotlarge = (delta_x_mean < QCcriteria.dxmax);
    
    % Get which pairs are not too far from cross-shore
    pwdissipation.data(i).lroughlyonx = (abs(dyondx) < QCcriteria.dyondxmax);
  
    %
    pwdissipation.data(i).bottomdepth_mean = mean(pwdissipation.data(i).bottomdepth, 1, 'omitnan');
    
    % Get pairs where mean bottom depth is sufficiently deep
    pwdissipation.data(i).lhmeandeepenough = ...
                (pwdissipation.data(i).bottomdepth_mean(1) > QCcriteria.hmin) & ...
                (pwdissipation.data(i).bottomdepth_mean(2) > QCcriteria.hmin);
    
end


%% Apply timeseries QC to flag data
%
% The single-parameter criteria are NOT used here to prevent
% doing computations of fe for certain pairs (because it is
% instructive to do calculations even when data are not appropriate
% to see what comes out of the calculation).

%
for i = 1:pwdissipation.Npairs
            
	%
    lflag_aux = ~pwdissipation.data(i).(fluxsource).lalltmsrsQC;
    
    %
    pwdissipation.data(i).Urms_avg(lflag_aux) = NaN;
    %
    pwdissipation.data(i).(fluxsource).Flux(lflag_aux, :) = NaN;
    pwdissipation.data(i).(fluxsource).dFluxdx(lflag_aux) = NaN;
    
end