function [r_ECI, v_ECI, extras] = PositionVelocityFromOrbitalElements(OE, mu)
% OrbitalElementsFromPositionVelocity Compute classical orbital elements
%   from position and velocity vectors in ECI frame.
%
% INPUTS:
%   OE         : struct with fields (units in brackets)
%                - a     : semi-major axis [m]
%                - e     : eccentricity [-]
%                - i     : inclination [rad]
%                - Omega : right ascension of ascending node, RAAN [rad]
%                - omega : argument of periapsis [rad]
%                - theta : true anomaly [rad]
%   mu (optional) : gravitational parameter [m^3/s^2]
%                   Default: getConstants().Earth.mu
%
% OUTPUTS:
%   r_ECI   : 3×n position vectors in ECI [m]
%   v_ECI   : 3×n velocity vectors in ECI [m/s]
%   extras  : struct with intermediate quantities:
%             - p        : semi-latus rectum [m]
%             - r        : radius magnitude [m]
%             - r_w      : position in perifocal (PQW) frame [m]
%             - v_w      : velocity in perifocal (PQW) frame [m/s]
%             - R_W2ECI  : 3×3×n DCMs from perifocal to ECI (ZXZ sequence)

if nargin<2
    mu      = getConstants().Earth.mu;  % Earth gravitational parameter
end

requiredFields = {'a','e','i','Omega','omega','theta'};
for k = 1:length(requiredFields)
    if ~isfield(OE, requiredFields{k})
        error('OE must contain field "%s"', requiredFields{k});
    end
end
n = length(OE.a);

extras.p = (OE.a*(1-OE.e.^2));
extras.r = extras.p./(1+OE.e.*cos(OE.theta));
extras.r_w =  [cos(OE.theta); sin(OE.theta);zeros(1,n)].*extras.r;
extras.v_w = [-sin(OE.theta); OE.e+cos(OE.theta);zeros(1,n)].*(mu./extras.p).^0.5;

extras.R_W2ECI = zeros(3,3,n);
r_ECI = zeros(3,n);
v_ECI = zeros(3,n);
for ii = 1:n
    extras.R_W2ECI(:,:,ii) = angle2dcm(-OE.omega(ii), -OE.i(ii), -OE.Omega(ii),'ZXZ');
    r_ECI(:,ii) = extras.R_W2ECI(:,:,ii)*extras.r_w(:,ii);
    v_ECI(:,ii) = extras.R_W2ECI(:,:,ii)*extras.v_w(:,ii);
end