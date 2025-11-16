function [rho, Az, El, rhoDot, AzDot, ElDot] = groundsStationOutput(GS,X, Phi_G)

r_ECI = X(1:3);
v_ECI = X(4:6);

lat = GS.LatDeg  * getConstants().Conversions.deg2rad;
lon_lst = GS.LonDeg  * getConstants().Conversions.deg2rad + Phi_G;
R_SEZ2ECI =  angle2dcm(0, -(pi/2 - lat), -lon_lst, 'XYZ');
r_h =  R_SEZ2ECI' * r_ECI;
rho_h = r_h -  [0;0;getConstants().Earth.Re];
rho = norm(rho_h);

El = asin(rho_h(3)/rho);
Az = atan2(rho_h(2), -rho_h(1));

rhoDot_ECI = v_ECI - cross([0;0;getConstants().Earth.omega], r_ECI);
rhoDot_h = R_SEZ2ECI' * rhoDot_ECI;
A = [-cos(Az)*cos(El),  rho*sin(Az)*cos(El),    rho*cos(Az)*sin(El);...
     sin(Az)*cos(El),   rho*cos(Az)*cos(El),    -rho*sin(Az)*sin(El);...
     sin(El),           0,                      rho*cos(El)];
rho_az_el_dot = A\rhoDot_h;
rhoDot = rho_az_el_dot(1);
AzDot = rho_az_el_dot(2);
ElDot = rho_az_el_dot(3);
