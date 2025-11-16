function Constants = getConstants()
%GETCONSTANTS  Returns a struct of useful constants for orbital mechanics.
%
%   OUTPUT:
%     Constants : struct containing conversion factors and celestial body data

%% Conversions
Constants.Conversions.Au2m   = 149597870700;
Constants.Conversions.m2Au   = 1/Constants.Conversions.Au2m;
Constants.Conversions.km2m   = 1e3;
Constants.Conversions.m2km   = 1e-3;
Constants.Conversions.mm2m   = 1e-3;
Constants.Conversions.years2s = 365.25*24*60*60;
Constants.Conversions.s2years = 1/Constants.Conversions.years2s;
Constants.Conversions.days2s = 24*60*60;
Constants.Conversions.s2days = 1/Constants.Conversions.days2s;
Constants.Conversions.hr2s   = 60*60;
Constants.Conversions.s2hrs  = 1/Constants.Conversions.hr2s;
Constants.Conversions.min2s  = 60;
Constants.Conversions.s2min  = 1/Constants.Conversions.min2s;
Constants.Conversions.deg2rad = pi/180;
Constants.Conversions.rad2deg = 180/pi;
Constants.Conversions.rev2rad = 2*pi;
Constants.Conversions.rad2rev = 1/Constants.Conversions.rev2rad;

%% Earth Constants
Constants.Earth.mu    = 398600 * Constants.Conversions.km2m^3; % m^3/s^2
Constants.Earth.omega = 7.292116e-5; % rad/s
Constants.Earth.Re = 6378 * Constants.Conversions.km2m; % m
Constants.Earth.J2 = 0.0010826; % Spherical harmonics terms J2
 
%% Sun Constants
Constants.Sun.mu = 132712000000 * Constants.Conversions.km2m^3; % m^3/s^2
Constants.Sun.Rs = 6.957e8;

%% Moon Constants 
Constants.Moon.mu = 4905 * Constants.Conversions.km2m^3; % m^3/s^2
Constants.Moon.Rm  = 1737 * Constants.Conversions.km2m;   % m
Constants.Moon.T  = 27.3 * Constants.Conversions.days2s; % s

%% Universal Constants
Constants.SpeedOfLight = 299792458;