function [pwdissipation] = pairwise_fe_calculate(pwdissipation)
%% [pwdissipation] = PAIRWISE_FE_CALCULATE(pwdissipation)
%
%   inputs
%       - pwdissipation: structure variable created by
%                        pairwise_fe_getdatapair.m
%
%   outputs
%       - pwdissipation: structure with additional fields relevant for the
%                        calculation of the friction factor (fe).
%
%
% PAIRWISE_FE_CALCULATE.m is a high-level function for computing the
% friction factor (fe) between an instrument pair in the ROXSI2022 dataset.
% The data from the pair is organized by the function
% pairwise_fe_getdatapair.m
%
%
% See also:
%   run_proc_pairwise_fe_allsites.m
%   pairwise_fe_getdatapair.m


%%
% --------------------------------------------
% --- COMPUTE TIMESERIES OF FRICTION FACTOR --
% --------------------------------------------

%
list_fields_aux = ["fromEcg", "fromFx"];

%
for i1 = 1:pwdissipation.Npairs
    %
    for i2 = 1:length(list_fields_aux)
        
        %
        pwdissipation.data(i1).(list_fields_aux(i2)).lposD = ...
               (-pwdissipation.data(i1).(list_fields_aux(i2)).dFluxdx) > 0;
        
        %
        pwdissipation.data(i1).(list_fields_aux(i2)).fe_timeseries = ...
                             -pwdissipation.data(i1).(list_fields_aux(i2)).dFluxdx ./ ...
                              (0.8*1025 .* pwdissipation.data(i1).Urms_avg.^3);        
                          
        % Compute time-mean
        pwdissipation.data(i1).(list_fields_aux(i2)).fe_mean = ...
            mean(pwdissipation.data(i1).(list_fields_aux(i2)).fe_timeseries, 'omitnan');
        
    end                
end


%%
% --------------------------------------------
% ------ COMPUTE PAIR-WISE CONSTANT/BULK -----
% ------------- FRICTION FACTORS -------------
% --------------------------------------------


%%

%
for i1 = 1:pwdissipation.Npairs
    %
    for i2 = 1:length(list_fields_aux)
        
        %
        var_1_aux = -pwdissipation.data(i1).(list_fields_aux(i2)).dFluxdx(:) ./ (0.8*1025);
        var_2_aux = pwdissipation.data(i1).Urms_avg(:).^3;
        
        %
        l_posD_aux = pwdissipation.data(i1).(list_fields_aux(i2)).lposD;

        %
        [pwdissipation.data(i1).(list_fields_aux(i2)).fctefit, ...
         pwdissipation.data(i1).(list_fields_aux(i2)).fcteerror] = ...
                                            regress(var_1_aux, ...
                                                    var_2_aux);
        
        % Correlation coefficient
        pwdissipation.data(i1).(list_fields_aux(i2)).rcoef = corr(var_1_aux, var_2_aux, 'rows', 'complete');
    
        
        % Correlation coefficient with just positive D
        if any(l_posD_aux)
            %
            [fe_bulk_aux, fe_bulk_error] = regress(var_1_aux(l_posD_aux), ...
                                                   var_2_aux(l_posD_aux));
            %
            rcoef_posD_aux = corr(var_1_aux(l_posD_aux), var_2_aux(l_posD_aux), 'rows', 'complete');
        else
            %
            fe_bulk_aux = NaN;
            fe_bulk_error = NaN;
            %
            rcoef_posD_aux = NaN;
        end
        %
        pwdissipation.data(i1).(list_fields_aux(i2)).posD.fctefit = fe_bulk_aux;
        pwdissipation.data(i1).(list_fields_aux(i2)).posD.fcteerror = fe_bulk_error;
        pwdissipation.data(i1).(list_fields_aux(i2)).posD.rcoef = rcoef_posD_aux;
        

    end
end