%% Define reference coordinates and direction
% for the different cartesian grids in ROXSI

%
clear


%% Coordinates and angles of reference for each array
%
% The angle is the arc (in clockwise degrees) from the
% geographical north to the offshore direction.

% ------------
% Asilomar
%
roxsigrid.Asilomar.latref = 36 +  (37/60) + (26.5187/3600);
roxsigrid.Asilomar.lonref = -121 -(56/60) - (25.19048/3600);
%
% roxsigrid.Asilomar.angleref = -23;    % angle of rotation (270 - 293)
roxsigrid.Asilomar.angleref = 293;

% ------------
% China Rock
%
roxsigrid.ChinaRock.latref = 36 +  (36/60) + (15.8927651999/3600);
roxsigrid.ChinaRock.lonref = -121 -(57/60) - (33.8133851999/3600);
%
% roxsigrid.ChinaRock.angleref = -15;    % angle of rotation (270 - 285)
roxsigrid.ChinaRock.angleref = 285;


%% Save roxsiGrid structure in the same folder as this script

%
fullnamepath = mfilename('fullpath');
indslash = strfind(fullnamepath, '/');

%
save(fullfile(fullnamepath(1:indslash(end)), 'ROXSI_xygrids.mat'), 'roxsigrid')

