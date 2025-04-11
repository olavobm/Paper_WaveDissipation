function spotterL2 = Spotter_mergeL2(fileA, varargin)
%% spotterL2 = SPOTTER_MERGEL2(fileA, varargin)
% 
%   inputs
%       - fileA: *.mat file with spotterL2 data structure.
%       - varargin: one (or many) *.mat files with spotterL2  structures.
% 
%   outputs
%       - spotterL2: a data structure with merged data from the input
%                    data files.
% 
% 
% SPOTTER_MERGEL2.m merges L2 Spotter data files from the ROXSI experiment,
% for the case when a Spotter had issues and was replaced by another.
% 
% SPOTTER_MERGEL2.m requires that L2 data structures from individual files
% only have data for the corresponding times when the Spotters were in the
% water. If they all have the same time grid (dtime), this function will
% fail.


%% Check the additional files

%
Nmore = length(varargin);
%
spotdata = cell(1, 1+Nmore);

%
if Nmore==0
    warning('Only one L2 file given in input. No merging performed')
    return
end


%% Put data from fileA in the first field of cell array

%
data_aux = load(fileA);
%
spotdata{1} = data_aux.spotterL2;


%% Load the rest of the data into the other fields of the cell

%
for i = 1:Nmore
    
    %
    data_aux = load(varargin{i});
    %
    spotdata{1+i} = data_aux.spotterL2;
    
end


%% Put spotdata in chronological order

%
dtime_lims = NaT(length(spotdata), 2);
dtime_lims.TimeZone = spotdata{1}.dtime.TimeZone;

%
for i = 1:length(spotdata)
    
    %
    dtime_lims(i, 1) = spotdata{i}.dtime(1);
    dtime_lims(i, 2) = spotdata{i}.dtime(end);
    
end

%
if any(dtime_lims(2:end, 1) == dtime_lims(1, 1))
    warning('dtime from the different files look the same.')
end

%
[~, ind_sort] = sort(dtime_lims(:, 1));
%
spotdata = spotdata(ind_sort);


%% Create the full time grid

%
dt_grid = diff(spotdata{1}.dtime(1:2));

%
dtime_grid = spotdata{1}.dtime(1) : dt_grid : spotdata{end}.dtime(end);
dtime_grid = dtime_grid(:);


%% Create output data structure by copying the structure from fileA

%
spotterL2 = spotdata{1};

%
spotterL2.dtime = dtime_grid;

%
list_fields = fieldnames(spotterL2);


%% Add/fill fields

% % %
% % list_refloc = ["X", "Y"];

%
for i1 = 1:length(spotdata)
    %
    ind_first = find(spotterL2.dtime == spotdata{i1}.dtime(1));
    ind_last = find(spotterL2.dtime == spotdata{i1}.dtime(end));
    %
    inds_fill = ind_first : 1 : ind_last;
    
    %
    if length(inds_fill)==1
        error('!!!')
    end
    
    %
    for i2 = 1:length(list_fields) 
        %
        if strcmp(list_fields{i2}, 'dtime')
            continue
        end
        %
        if strcmp(list_fields{i2}, 'SN')
            if i1>1
                spotterL2.SN = strcat(spotterL2.SN, ", ", spotdata{i1}.SN);
            end
        end

        %
        if size(spotdata{i1}.(list_fields{i2}), 1)==length(spotdata{i1}.dtime)
            %
            if i1==1
                spotterL2.(list_fields{i2}) = NaN(length(spotterL2.dtime), ...
                                                   size(spotdata{i1}.(list_fields{i2}), 2));
            end
            %
            spotterL2.(list_fields{i2})(inds_fill, :) = spotdata{i1}.(list_fields{i2});
        end
    end
    
end


%% Similar as above, but for subfields of location

%
list_fields = fieldnames(spotterL2.location);

%
for i1 = 1:length(spotdata)
    %
    ind_first = find(spotterL2.dtime == spotdata{i1}.dtime(1));
    ind_last = find(spotterL2.dtime == spotdata{i1}.dtime(end));
    %
    inds_fill = ind_first : 1 : ind_last;
    
    %
    for i2 = 1:length(list_fields) 
        %
        if strcmp(list_fields{i2}, 'dtime')
            continue
        end
        
        %
        if size(spotdata{i1}.location.(list_fields{i2}), 1)==length(spotdata{i1}.dtime)
            %
            if i1==1
                spotterL2.location.(list_fields{i2}) = NaN(length(spotterL2.dtime), ...
                                                           size(spotdata{i1}.location.(list_fields{i2}), 2));
            end
            %
            spotterL2.location.(list_fields{i2})(inds_fill, :) = spotdata{i1}.location.(list_fields{i2});
        end
    end
    
end




