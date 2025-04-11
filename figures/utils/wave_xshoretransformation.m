function modeloutput = wave_xshoretransformation(x, h, dx, freq0, Hrms0, angle0, gamma, B)
%% modeloutput = WAVE_XSHORETRANSFORMATION(x, h, dx, freq0, Hrms0, angle0, gamma, B)
%
%   inputs
%       - x: cross-shore location (in meters).
%       - h: bottom depth (positive, in meters).
%       - dx: grid resolution (in meters). If given and different than x,
%             then linearly interpolate bottom depth.
%       - freq0: wave frequency (in Hz).
%       - Hrms0: RMS wave height (in meters) at the offshore boundary h(1).
%       - angle0: angle of incidence (in radians) at the offshore boundary.
%       - gamma: depth-limited wave breaking coefficient.
%       - B: O(1) breaker coefficient.
%
%   outputs
%       - modeloutput: structure variable with all data and parameters in
%                      the wave transformation model.
%
%
% WAVE_XSHORETRANSFORMATION.m computes the cross-shore surface gravity wave
% transformation as in Thornton and Guza (1983), for 0 bottom friction.


%% Define parameters if not given as inputs

%
if ~exist('angle0', 'var') || isempty(angle0)
    angle0 = 0;
end

%
if ~exist('gamma', 'var') || isempty(gamma)
    gamma = 0.45;
end

%
if ~exist('B', 'var') || isempty(B)
    B = 1;
end


%% Check grid resolution

%
if isempty(dx)
    dx = diff(x(1:2));
end

%
if dx == diff(x(1:2))
    
    %
    xinterp = x(:);
    hinterp = h(:);
%
else
    %
    xinterp = x(1) : dx : x(end);
    xinterp = xinterp(:);
    
    %
    hinterp = interp1(x, h, xinterp);
end


%% Pre-allocate space for output variables

%
modeloutput.parameters.frequency = freq0;
modeloutput.parameters.gamma = gamma;
modeloutput.parameters.B = B;

%
N = length(hinterp);

%
modeloutput.dx = dx;
modeloutput.x = xinterp(:);
modeloutput.h = hinterp(:);

%
modeloutput.Hrms = NaN(N, 1);
%
modeloutput.epsilon_b = modeloutput.Hrms;

%
list_othervars = ["k", "cp", "angleinc", "cg", "Energy", "Ecgx"];
%
for i = 1:length(list_othervars)
    modeloutput.(list_othervars(i)) = modeloutput.Hrms;
end


%% See if h < h_cutoff near onshore boundary, and remove from
% the integration (h_cutoff may be greater than 0 to avoid
% errors when computing wavenumber at very shallow water)

%
h_cutoff = 0;

%
lhtooshallow = (hinterp <= h_cutoff);
lhnan = isnan(hinterp);
%
lhbad = (lhtooshallow | lhnan);

%
ind_hbad = find(lhbad, 1, 'first');
%
if isempty(ind_hbad)    % if there is no bad bottom depth
    ind_hbad = length(hinterp) + 1;
end

%
ind_hok = 1:(ind_hbad-1);


%% Cross-shore integrate energy density equation (Thornton and Guza, 1983).

%
[Hrmsout, ...
 epsilonout, ...
 paramsout] = wave_transformation(dx, h, ...
                                  freq0, Hrms0, angle0, ...
                                  gamma, B);


                             
%% Now check for complex Hrmsout (when dissipation is
% larger than the incoming cross-shore energy flux)
    
%
lrealH = (Hrmsout == real(Hrmsout));
lcomplexH = ~lrealH;

%
if any(lcomplexH)
    %
    Hrmsout(lcomplexH) = 0;    % this removes complex Attribute of Hrmsout
    epsilonout(lcomplexH, :) = 0;
end


%% Organize output variables

%
modeloutput.Hrms(ind_hok) = Hrmsout;
%
modeloutput.epsilon_b(ind_hok) = epsilonout(:, 1);

%
modeloutput.k(ind_hok) = paramsout.k(:);
modeloutput.cp(ind_hok) = paramsout.cp(:);
modeloutput.angleinc(ind_hok) = paramsout.angle(:);
modeloutput.cg(ind_hok) = paramsout.cg(:);
modeloutput.Energy(ind_hok) = paramsout.E(:);
modeloutput.Ecgx(ind_hok) = paramsout.Ecgx(:);


%%
% ------------------------------------------------------
% ----------------- EMBEDDED FUNCTIONS -----------------
% ------------------------------------------------------

%
function x = rho()
    x = 1025;
end
%
function x = g()
    x = 9.81;
end


%%

function [Hrmsout, epsilonout, paramsout] = wave_transformation(dx, h, freq0, Hrms0, angle0, gamma, B)
%%
%
% WAVE_TRANSFORMATION_ROUGHNESS.m integrates the energy density equation in
% the cross-shore direction similar to Thornton and Guza (1983), a.k.a.
% TG83.

    %% Get number of grid points
    
    %
    Ngpt = length(h);
    
    
    %% Compute wavenumber and group velocity at all grid points
    
    %
    k = wave_freqtok(freq0, h);
    
    %
    cg = wave_cg(k(:), h(:));
    
    %
    cp = 2*pi*freq0 ./ k;
    
    
    %% Compute angle of incidence from Snell's law (note
    % the angle at all locations can be computed from
    % the angle at the boundary)
    
    %
    sin_theta = (cp./cp(1)) * sin(angle0);
    
    %
    theta = asin(sin_theta);
    
    %
    cos_theta = cos(theta);
    
    
    %% Cross-shore energy flux density at the offshore boundary
    
    %
    Ecgx0 = (1/8) * rho * g * (Hrms0^2) * cg(1) * cos_theta(1);
    
    
    %% Dissipation by wave breaking at the offshore boundary
    
    %
    epsi0_breaking = dissipation_breaking(Hrms0, freq0, h(1), gamma, B);
                 
    
    %% Pre-allocate arrays for variables computed in the loop below
    
    %
    nan_array = NaN(Ngpt, 1);
    %
    Hrmsout = nan_array;
    Energy_all = nan_array;
    xshoreFlux_all = nan_array;
    dissipation_breaking_all = nan_array;
    
    %
    Hrmsout(1) = Hrms0;
    %
    Energy_all(1) = (1/8) * rho * g * (Hrms0^2);
    %
    xshoreFlux_all(1) = Ecgx0;
    %
    dissipation_breaking_all(1) = epsi0_breaking;
    
    
    %% Loop over grid points and cross-shore integrate energy equation
    
    % Values at the boundary used in the loop below
    Ecgx_old = Ecgx0;
    epsilon_breaking_old = epsi0_breaking;
    
    %
    for igp = 2:Ngpt
        
        % Compute new E*cg*cos(theta) -- Equation (43) in TG83
        Ecxgnew = Ecgx_old - dx*epsilon_breaking_old;
    
        % Compute new Hrms
        Energy_new = (Ecxgnew./(cg(igp)*cos_theta(igp)));
        Hrms_new = sqrt(Energy_new * (8/(rho*g)));
        
        % Update variables for next loop iteration
        Ecgx_old = Ecxgnew;
        epsilon_breaking_old = dissipation_breaking(Hrms_new, freq0, h(igp), gamma, B);
     
        % Save variables in output arrays
        Hrmsout(igp) = Hrms_new;
        xshoreFlux_all(igp) = Ecxgnew;
        Energy_all(igp) = Energy_new;
        dissipation_breaking_all(igp) = epsilon_breaking_old;
    end
    
    
    %% Define additional output variables
    
    %
    epsilonout = dissipation_breaking_all;
    
    %
    paramsout.gamma = gamma;
    paramsout.B = B;
    %
    paramsout.h = h(:);
    %
    paramsout.k = k;
    paramsout.angle = theta;
    paramsout.cg = cg;
    paramsout.cp = cp;
    %
    paramsout.E = Energy_all;
    paramsout.Ecgx = xshoreFlux_all;


end


% Dissipation by wave breaking -- Equation (26) in TG83
function epsilon = dissipation_breaking(Hrms, freq, h, gamma, B)
    
    %
    epsilon = (3*sqrt(pi)/16) * (rho*g) * ...
              ((B^3)*freq) * ...
              (Hrms^5)/((gamma^2)*(h^3)) * ...
              (1 - (1 + (Hrms/(gamma*h))^2).^(-5/2) );

end



%
end
