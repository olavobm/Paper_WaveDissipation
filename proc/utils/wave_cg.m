function [cg, kH] = wave_cg(k, H)
%% [cg, kH] = WAVE_CG(k, H)
%
%   inputs
%       - k: wavenumber (in radians per meter).
%       - H: bottom depth (in meters, positive).
%
%   outputs
%       - cg: group velocity (in m/s).
%
%
% WAVE_CG.m computes the group velocity of surface gravity waves
% using linear wave theory. Requires wave_cp.m
%
% See also:
%   wave_cp.m
%
% Olavo Badaro Marques


%
kH = k .* H;

%
cp = wave_cp(k, H);

%
cg = cp .* 0.5 .* (1 + (2*kH./sinh(2*kH)));

