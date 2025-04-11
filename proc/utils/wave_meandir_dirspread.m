function [theta1, sigma1, theta2, sigma2] = wave_meandir_dirspread(a1, b1, a2, b2)
%% [theta1, sigma1, theta2, sigma2] = WAVE_MEANDIR_DIRSPREAD(a1, b1, a2, b2)
%
%   inputs
%       - a1, a2, b1, b2: directional moments.
%
%   outputs
%       - theta1: mean direction based on a1 and b1.
%       - sigma1: directional spread based on a1 and b1.
%       - theta2: mean direction based on a2 and b2.
%       - sigma2: directional spread based on a1, b1, a2, and b2.
%
%
% WAVE_MEANDIR_DIRSPREAD.m computes meand direction and directional spread
% from the first and second wave directional moments. All results have
% units of radians.
%
% Along with the definitions in the lower level functions computing
% directional moments (i.e. See also below), theta1 gives the mean
% direction (in radians) of where waves propagate to.
%
% theta2 is oblivious about whether waves are propagating
% to or from the theta2 direction.
%
%
% See also:
%   wave_dirmoments_xyz.m
%   wave_dirmoments_uvp.m

%% Compute first set of mean direction and directional spread

%
theta1 = atan2(b1, a1);

% Herbers et al. (1999)
sigma1 = sqrt(2 * (1 - a1.*cos(theta1) - b1.*sin(theta1)));

% % Herbers et al. (2012) -- apart from machine precision, this is the same as sigma1 
% sigma1 = sqrt(2 * (1 - sqrt(a1.^2 + b1.^2)));


%% Compute second set of mean direction and directional spread

%
theta2 = 0.5*atan2(b2, a2);

% Herbers et al. (1999)
sigma2 = sqrt(0.5 * (1 - a2.*cos(2*theta1) - b2.*sin(2*theta1)));

% % Herbers et al. (2012) -- this is NOT the same as sigma2
% sigma2 = sqrt(0.5 * (1 - sqrt(a2.^2 + b2.^2)));


