function dataL2 = roxsi_retimeL2(dataL2in, dtime_lims)
%% dataL2 = ROXSI_RETIMEL2(dataL2in, dtime_lims)
% 
%   inputs
%       - dataL2in: L2 data structure from ROXSI.
%       - dtime_lims: datetime limits to grab the data.
% 
%   outputs
%       - dataL2: a data structure with merged data from the input
%                    data files.
% 
% 
% ROXSI_RETIMEL2.m grabs L2 data within dtime_lims. This function is 
% based on Spotter_retimeL2.m, because it should work for any L2 data
% structure in ROXSI.
%
% As it is, ROXSI_RETIMEL2.m only works with variables that are 1D or
% 2D arrays. For 3D arrays (e.g. directional spectra), something more
% complicated is required in the code.
%
% See also:
%   Spotter_retimeL2.m



%% Create time grid of the output

%
if isempty(dtime_lims.TimeZone)
    dtime_lims.TimeZone = "America/Los_Angeles";
end

%
dt_grid = diff(dataL2in.dtime(1:2));

%
dtime_grid = dtime_lims(1) : dt_grid : dtime_lims(2);



%% Create output data structure by copying the structure from fileA

%
dataL2 = dataL2in;
%
dataL2.dtime = dtime_grid;


%%

%
ind_get_first = find(dataL2in.dtime >= dataL2.dtime(1), 1, 'first');
ind_get_last = find(dataL2in.dtime <= dataL2.dtime(end), 1, 'last');

%
if isempty(ind_get_first)
    %
    warning(['No data is available at the given time grid points. ' ...
             'Output structure will only have NaN for the data.'])
end



%%

%
ind_fill_first = find(dataL2.dtime >= dataL2in.dtime(ind_get_first), 1, 'first');
ind_fll_last = find(dataL2.dtime <= dataL2in.dtime(ind_get_last), 1, 'last');


%%

%
inds_fill = ind_fill_first : 1 : ind_fll_last;
%
inds_get = ind_get_first : 1 : ind_get_last;


%% Add/fill fields of the outer structure

%
list_fields = fieldnames(dataL2in);

%
for i = 1:length(list_fields)
    
    %
    if strcmp(list_fields{i}, 'dtime')
        continue
    end

    %
    if size(dataL2in.(list_fields{i}), 1) == length(dataL2in.dtime)
        %
        dataL2.(list_fields{i}) = NaN(length(dataL2.dtime), ...
                                         size(dataL2in.(list_fields{i}), 2));
                                     
        %
        dataL2.(list_fields{i})(inds_fill, :) = dataL2in.(list_fields{i})(inds_get, :);
    end
end


%% Similar as above, but for subfields of location (in Spotter data)

%
if isfield(dataL2in, 'location')

    %
    list_fields = fieldnames(dataL2in.location);

    %
    for i = 1:length(list_fields) 
        %
        if strcmp(list_fields{i}, 'dtime')
            continue
        end

        %
        if size(dataL2in.location.(list_fields{i}), 1)==length(dataL2in.dtime)
            %
            dataL2.location.(list_fields{i}) = NaN(length(dataL2.dtime), ...
                                                   size(dataL2in.location.(list_fields{i}), 2));
            %
            dataL2.location.(list_fields{i})(inds_fill, :) = dataL2in.location.(list_fields{i})(inds_get, :);
        end
    end

end


