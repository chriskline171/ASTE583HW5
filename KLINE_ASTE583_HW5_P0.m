clear all
%% Ground station definitions
Stations.One.LatDeg    = -35.398333;
Stations.One.LonDeg    = 148.981944;

Stations.Two.LatDeg    = 40.427222;
Stations.Two.LonDeg    = 355.749444;

Stations.Three.LatDeg  = 35.247164;
Stations.Three.LonDeg  = 243.205;

%% Load and preprocess raw data
rawData = load('RawMeasurements.txt');      % [t, stnID, rho_km, rhodot_kmps]
rawData(:,3:4) = rawData(:,3:4)*getConstants().Conversions.km2m; % convert to m/m/s

%% Initial OE guess
OE0.a     = 10000*getConstants().Conversions.km2m;
OE0.e     = 0.001;
OE0.i     = 40*getConstants().Conversions.deg2rad;
OE0.Omega = 80*getConstants().Conversions.deg2rad;
OE0.omega = 40*getConstants().Conversions.deg2rad;
OE0.theta = 0;
[r0_ECI, v0_ECI] = PositionVelocityFromOrbitalElements(OE0);

%% Build augmented initial state (r,v,STM)
X0 = [r0_ECI; v0_ECI; reshape(eye(6),36,1)];

t0 = 0;
tf = 24*getConstants().Conversions.hr2s;

% force ode45 output at measurement times
tspan = [t0, rawData((rawData(:,1)<tf & rawData(:,1)>t0),1)', tf];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
phi_G0 = 170*getConstants().Conversions.deg2rad;

%% Batch filter (Gauss–Newton on initial state)
delta_x = 0;            % correction vector
firtIteration = true;
kk = 0;

while max(abs(delta_x))>1e-3 || firtIteration
    firtIteration = false;
    X0(1:6) = X0(1:6) + delta_x;      % update initial state guess

    [t,X] = ode45(@twobodyEOMJ2WithSTM, tspan, X0, options);

    % determine if endpoints contain measurements
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

    delta_y = nan(length(indxStart:indxEnd),4);  % [t, stnID, drho, drhodot]
    phi_G = getConstants().Earth.omega*t + phi_G0;
    stationNames = fieldnames(Stations);

    nn = 1;
    rawDataStartIndx = find(rawData(:,1)>=t0,1,'first')-1;

    HTransposeH = zeros(6,6);
    HTransposeDeltaY = zeros(6,1);

    %% Loop over measurement points
    for ii = indxStart:indxEnd
        delta_y(nn,1) = t(ii);
        stnID = rawData(rawDataStartIndx+nn,2);
        delta_y(nn,2) = stnID;

        % modelled measurement
        [rho,~,~,rhoDot,~,~] = ...
            groundsStationOutput(Stations.(stationNames{stnID}), ...
                                 X(ii,:)', phi_G(ii));

        % residuals
        delta_y(nn,3) = rawData(rawDataStartIndx+nn,3) - rho;
        delta_y(nn,4) = rawData(rawDataStartIndx+nn,4) - rhoDot;

        % measurement sensitivity wrt current (r,v)
        Htld = calculateH_rho_rhoDot_wrt_r_v( ...
                 X(ii,:)', Stations.(stationNames{stnID}), phi_G(ii));

        % STM at this time
        Phi = reshape(X(ii,7:end),6,6);
        H = Htld * Phi;    % sensitivity wrt initial state

        % accumulate normal eqns
        HTransposeH      = HTransposeH      + H'*H;
        HTransposeDeltaY = HTransposeDeltaY + H'*delta_y(nn,3:4)';

        nn = nn + 1;
    end

    delta_x = HTransposeH \ HTransposeDeltaY;   % Gauss-Newton update
    kk = kk + 1;
end

%% Final post-fit difference from OE-based initial state
delta_x_postfit = X0(1:6) - [r0_ECI; v0_ECI];
