function ind_match = find_ind_matchstring(strarray, list_str)
%% ind_match = FIND_IND_MATCHSTRING(strarray, list_str)
%
%   inputs
%       - strarray: a vector of strings.
%       - list_str: an array of strings to search in strarray.
%
%   outputs
%       - ind_match: indices of strarray corresponding
%                    to strings in list_str.
%
%
% FIND_IND_MATCHSTRING.m search in strarray for the strings specified
% in list_str. Return the corresponding indices.


%% Check inputs

if ~isvector(strarray)
    error('First input must be a vector')
end


%% Get size of list_str so that output has the same size and
% pre-allocate a column vector for an intermediate variable

%
size_input = size(list_str);

%
vec_ind_match = NaN([numel(list_str), 1]);


%% Loop over list_str and find the indices of
% the corresponding strings in strarray

%
vec_list_str = list_str(:);

%
for i = 1:length(vec_list_str)
    
    %
    ind_match_aux = find(strcmp(strarray, vec_list_str(i)));
    
    %
    if isscalar(ind_match_aux)
        vec_ind_match(i) = ind_match_aux;
    end
    
end

%
ind_match = reshape(vec_ind_match, size_input);
