function cp = wave_cp(k, H)
%% cg = WAVE_CP(k, H)
%
%   inputs
%       - k: wavenumber (in radians per meter).
%       - H: bottom depth (in meters, positive).
%
%   outputs
%       - cp: phase speed (in m/s).
%
%
% WAVE_CP.m computes phase speed of surface gravity waves
% using linear wave theory.
%
%
% Olavo Badaro Marques.


%
g = 9.8;

%
cp = sqrt(g * tanh(k.*H) ./ k);