%% Ground station definitions
Stations.One.LatDeg    = -35.398333;
Stations.One.LonDeg    = 148.981944;

Stations.Two.LatDeg    = 40.427222;
Stations.Two.LonDeg    = 355.749444;

Stations.Three.LatDeg  = 35.247164;
Stations.Three.LonDeg  = 243.205;

%% Load and preprocess raw data
rawData = load('TruthMeasurements.txt');      % [t, stnID, rho_km, rhodot_kmps]
rawData(:,3:4) = rawData(:,3:4)*getConstants().Conversions.km2m; % convert to m/m/s

%% Initial OE 
OE0.a     = 10000*getConstants().Conversions.km2m;
OE0.e     = 0.001;
OE0.i     = 40*getConstants().Conversions.deg2rad;
OE0.Omega = 80*getConstants().Conversions.deg2rad;
OE0.omega = 40*getConstants().Conversions.deg2rad;
OE0.theta = 0;
[r0_ECI, v0_ECI] = PositionVelocityFromOrbitalElements(OE0);

%% A priori data
sigma_r = 1*getConstants().Conversions.km2m;
sigma_v = 1;
P0 = [sigma_r^2*eye(3), zeros(3);
     zeros(3), sigma_v^2*eye(3)];
sigma_rho = 1;
sigma_rhoDot = 1*getConstants().Conversions.mm2m;
R = [sigma_rho^2, 0;
    0, sigma_rhoDot^2];
%% Build augmented initial state (r,v,STM)
X0 = [r0_ECI; v0_ECI; reshape(eye(6),36,1)];

t0 = 0;
tf = 24*getConstants().Conversions.hr2s;

% force ode45 output at measurement times
tspan = [t0, rawData((rawData(:,1)<tf & rawData(:,1)>t0),1)', tf];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
phi_G0 = 170*getConstants().Conversions.deg2rad;

%% Batch filter (Gauss–Newton on initial state)
delta_x_hat = zeros(6,1);            % correction vector
delta_x_bar = delta_x_hat;
firtIteration = true;
kk = 0;
maxIterLimit = 10;
while delta_x_hat'*(P0\delta_x_hat)>1e-15 || firtIteration
    firtIteration = false;
    X0(1:6) = X0(1:6) + delta_x_hat;      % update initial state guess

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

    InfoMat = inv(P0);
    InfoVec = P0\delta_x_bar;

    %% Loop over measurement points
    for ii = indxStart:indxEnd
        delta_y(nn,1) = t(ii);
        stnID = rawData(rawDataStartIndx+nn,2);
        delta_y(nn,2) = stnID;

        % modelled measurement
        [rho,~,~,rhoDot,~,~] = groundsStationOutput(Stations.(stationNames{stnID}),X(ii,:)', phi_G(ii));

        % residuals
        delta_y(nn,3) = rawData(rawDataStartIndx+nn,3) - rho;
        delta_y(nn,4) = rawData(rawDataStartIndx+nn,4) - rhoDot;

        % measurement sensitivity wrt current (r,v)
        Htld = calculateH_rho_rhoDot_wrt_r_v(X(ii,:)', Stations.(stationNames{stnID}), phi_G(ii));

        % STM at this time
        Phi = reshape(X(ii,7:end),6,6);
        H = Htld * Phi;    % sensitivity wrt initial state
        
        % accumulate normal eqns
        InfoMat = InfoMat + H'*(R\H);
        InfoVec = InfoVec + H'*(R\delta_y(nn,3:4)');

        nn = nn + 1;
    end

    delta_x_hat = InfoMat \ InfoVec;   % Gauss-Newton update
    delta_x_bar = delta_x_bar - delta_x_hat; % Update a priori state deviation
    kk = kk + 1;
    if kk>maxIterLimit
         warning(['Batch filter did not converge within %d iterations.\n' ...
             'Last step Mahalanobis norm: %.3e\n'], ...
             maxIterLimit, delta_x_hat'*(P0\delta_x_hat));
        break
    end
end
kk
%% Final post-fit difference from OE-based initial state
delta_x_postfit = X0(1:6) - [r0_ECI; v0_ECI];
