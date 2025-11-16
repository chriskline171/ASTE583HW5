function [OE, extras] = OrbitalElementsFromPositionVelocity(r_ECI, v_ECI, mu)
% OrbitalElementsFromPositionVelocity Compute classical orbital elements
%   from position and velocity vectors in ECI frame.
%
% INPUTS:
%   r_ECI : 3xn vector, position in ECI frame [m]
%   v_ECI : 3xn vector, velocity in ECI frame [m/s]
%   mu    : gravitational parameter [m^3/s^2]
%
% OUTPUTS:
%   OE     : structure containing classical orbital elements
%       a       : semi-major axis [m]
%       e       : eccentricity [unitless]
%       i       : inclination [rad]
%       Omega   : right ascension of ascending node (RAAN) [rad]
%       omega   : argument of periapsis [rad]
%       theta   : true anomaly [rad]
%   extras : structure containing useful intermediate values
%       h_ECI   : specific angular momentum vector [m^2/s]
%       e_ECI   : eccentricity vector [unitless]
%       p       : semi-latus rectum [m]
%       N_ECI   : node vector [unitless]
%       vr      : radial velocity [m/s]

if nargin<3
    mu      = getConstants().Earth.mu;  % Earth gravitational parameter
end
if (size(r_ECI,1) ~=3 || size(v_ECI,1) ~=3) || (size(r_ECI,2) ~= size(v_ECI,2))
    error('Length of input vectors must be 3-by-n')
else
    n = size(r_ECI,2);
end


% Compute angular momentum vector
extras.h_ECI = cross(r_ECI, v_ECI,1);

% Compute eccentricity vector
extras.e_ECI = cross(v_ECI, extras.h_ECI, 1)/mu - r_ECI./vecnorm(r_ECI,2, 1);
OE.e = vecnorm(extras.e_ECI, 2, 1);

% Compute semi-latus rectum and semi-major axis
extras.p = vecnorm(extras.h_ECI, 2, 1).^2/mu;
OE.a = extras.p./(1-OE.e.^2);

% Compute inclination
OE.i = acos(extras.h_ECI(3,:)./vecnorm(extras.h_ECI, 2, 1));


% Compute ascending node vector
extras.N_ECI = cross([zeros(2,n); ones(1,n)], extras.h_ECI,1);

% Compute right ascension of ascending node (RAAN)
N_ECI_Y_Positive_indx = extras.N_ECI(2,:)>=0;
OmegaNoCorrection = acos(extras.N_ECI(1,:)./vecnorm(extras.N_ECI(:,:),2,1));
OE.Omega(N_ECI_Y_Positive_indx) = OmegaNoCorrection(N_ECI_Y_Positive_indx);
OE.Omega(~N_ECI_Y_Positive_indx) = 2*pi-OmegaNoCorrection(~N_ECI_Y_Positive_indx);

% Compute argument of periapsis
e_ECI_Z_Positive_indx = extras.e_ECI(3,:)>=0;
omegaNoCorrection = acos(dot(extras.N_ECI(:,:),extras.e_ECI(:,:),1)./(vecnorm(extras.N_ECI(:,:),2,1).*vecnorm(extras.e_ECI(:,:),2,1)));
OE.omega(e_ECI_Z_Positive_indx) = omegaNoCorrection(e_ECI_Z_Positive_indx);
OE.omega(~e_ECI_Z_Positive_indx) = 2*pi-omegaNoCorrection(~e_ECI_Z_Positive_indx);


% Compute radial velocity
extras.vr = dot(r_ECI, v_ECI,1)./vecnorm(r_ECI,2,1);

% Compute true anomaly
vr_Positive_indx = extras.vr >= 0;
thetaNoCorrection =  acos(dot(extras.e_ECI, r_ECI,1)./(vecnorm(extras.e_ECI,2,1).*vecnorm(r_ECI,2,1)));
OE.theta(vr_Positive_indx) = thetaNoCorrection(vr_Positive_indx);
OE.theta(~vr_Positive_indx) = 2*pi - thetaNoCorrection(~vr_Positive_indx);


elipticalOrCircularOrbit = OE.e<=1;
extras.E(elipticalOrCircularOrbit) = atan2((1 - OE.e(elipticalOrCircularOrbit)).^0.5 .* sin(OE.theta(elipticalOrCircularOrbit)/2), (1 + OE.e(elipticalOrCircularOrbit)).^0.5 .* cos(OE.theta(elipticalOrCircularOrbit)/2) ) * 2;
extras.E(~elipticalOrCircularOrbit) = NaN;
% Mean Anomaly
extras.M(elipticalOrCircularOrbit) = extras.E(elipticalOrCircularOrbit) - OE.e(elipticalOrCircularOrbit).*sin(extras.E(elipticalOrCircularOrbit));
extras.M(~elipticalOrCircularOrbit) = NaN;
% Mean motion
extras.n = (mu./(OE.a.^3)).^0.5;