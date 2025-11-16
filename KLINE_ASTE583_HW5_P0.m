Stations.One.LatDeg = -35.398333;
Stations.One.LonDeg = 148.981944;

Stations.Two.LatDeg = 40.427222;
Stations.Two.LonDeg = 355.749444;

Stations.Three.LatDeg = 35.247164;
Stations.Three.LonDeg = 243.205;

% Test filter with Raw measurements
rawData =  load('RawMeasurements.txt');
rawData(:,3:4) = rawData(:,3:4)*getConstants().Conversions.km2m;
OE0.a = 10000*getConstants().Conversions.km2m;
OE0.e = 0.001;
OE0.i = 40 *getConstants().Conversions.deg2rad;
OE0.Omega = 80 *getConstants().Conversions.deg2rad;
OE0.omega = 40 *getConstants().Conversions.deg2rad;
OE0.theta = 0;
[r0_ECI, v0_ECI] = PositionVelocityFromOrbitalElements(OE0);


X0 = [r0_ECI; v0_ECI; reshape(eye(6,6),6*6,1)] ;
t0 = 0;
tf = 24*getConstants().Conversions.hr2s;
tspan = [t0, rawData((rawData(:,1)<tf & rawData(:,1)>t0),1)', tf];
% Integration options
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10) ;
phi_G0 = 170 * getConstants().Conversions.deg2rad;

%% Filter Algorthm
delta_x = 0; 
firtIteration = true;
kk = 0;
while max(abs(delta_x))>1e-3 || firtIteration == true 
    firtIteration = false;
    X0(1:6) = X0(1:6)+delta_x;
    [t,X] = ode45(@twobodyEOMJ2WithSTM, tspan, X0, options);
    if any(t0 == rawData(:,1))
        indxStart = 1;
    else
        indxStart = 2;
    end
    if any(tf == rawData(:,1))
        indxEnd = length(t);
    else
        indxEnd = length(t)-1;
    end
    delta_y = nan(length(indxStart:indxEnd), 4);
    phi_G = getConstants().Earth.omega*t+phi_G0;
    stationNames = fieldnames(Stations);
    nn = 1;
    rawDataStartIndx = find(rawData(:,1)>=t0, 1,'first')-1;
    HTransposeH = zeros(6,6);
    HTransposeDeltaY = zeros(6,1);
    for ii = indxStart:indxEnd
        delta_y(nn,1) = t(ii);
        delta_y(nn,2) = rawData(rawDataStartIndx+nn,2);
        [rho, ~, ~, rhoDot, ~, ~] = groundsStationOutput(Stations.(stationNames{rawData(rawDataStartIndx+nn,2)}), X(ii,:)', phi_G(ii));
        delta_y(nn,3) = rawData(rawDataStartIndx+nn,3) - rho;
        delta_y(nn,4) = rawData(rawDataStartIndx+nn,4) - rhoDot;
        Htld  = calculateH_rho_rhoDot_wrt_r_v(X(ii,:)', Stations.(stationNames{rawData(rawDataStartIndx+nn,2)}), phi_G(ii));
        H = Htld* reshape(X(ii,7:end), 6,6);
        HTransposeH = HTransposeH + H'*H;
        HTransposeDeltaY = HTransposeDeltaY + H'*delta_y(nn,3:4)';
        nn = nn+1;
    end
    
    delta_x = HTransposeH\HTransposeDeltaY;
    kk = kk+1;
end

delta_x_postfit = X0(1:6) - [r0_ECI; v0_ECI]