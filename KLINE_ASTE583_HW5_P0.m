clc; clear all; close all

%% Ground station definitions
Stations.One.LatDeg    = -35.398333;
Stations.One.LonDeg    = 148.981944;

Stations.Two.LatDeg    = 40.427222;
Stations.Two.LonDeg    = 355.749444;

Stations.Three.LatDeg  = 35.247164;
Stations.Three.LonDeg  = 243.205;

%% Load and preprocess raw data
measData = load('TruthMeasurements.txt');      % [t, stnID, rho_km, rhodot_kmps]
measData(:,3:4) = measData(:,3:4)*getConstants().Conversions.km2m; % convert to m/m/s

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
tspan = [t0, measData((measData(:,1)<tf & measData(:,1)>t0),1)', tf];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
phi_G0 = 170*getConstants().Conversions.deg2rad;

%% Batch filter (Gauss–Newton on initial state)
delta_x_hat = [0.541708589833822;0.215276743995227;-0.0246521235406893;-5.6364922705741e-5;-0.000394199357037521;0.00051130212438399]*getConstants().Conversions.km2m;            % correction vector
delta_x_bar = zeros(size(delta_x_hat));
firstIteration = true;
kk = 0;
maxIterLimit = 15;
MalhonbisDist =  delta_x_hat'*(P0\delta_x_hat);
while MalhonbisDist>1e-10 || firstIteration
    firstIteration = false;
    kk = kk + 1;
    X0(1:6) = X0(1:6) + delta_x_hat;      % update initial state guess
    IterationData{kk}.X0 = X0(1:6);
    IterationData{kk}.delta_x_hat = delta_x_hat;
    IterationData{kk}.MalhonbisDist = MalhonbisDist;
    IterationData{kk}.delta_x_bar = delta_x_bar;

    [t,X] = ode45(@twobodyEOMJ2WithSTM, tspan, X0, options);

    % determine if endpoints contain measurements
    if any(t0 == measData(:,1))
        indxStart = 1;
    else
        indxStart = 2;
    end
    if any(tf == measData(:,1))
        indxEnd = length(t);
    else
        indxEnd = length(t)-1;
    end

    delta_y = nan(length(indxStart:indxEnd),4);  % [t, stnID, drho, drhodot]
    epsilon_hat= nan(length(indxStart:indxEnd),2);
    phi_G = getConstants().Earth.omega*t + phi_G0;
    stationNames = fieldnames(Stations);

    nn = 1;
    rawDataStartIndx = find(measData(:,1)>=t0,1,'first')-1;

    InfoMat = inv(P0);
    InfoVec = P0\delta_x_bar;

    %% Loop over measurement points
    for ii = indxStart:indxEnd
        delta_y(nn,1) = t(ii);
        stnID = measData(rawDataStartIndx+nn,2);
        delta_y(nn,2) = stnID;

        % modelled measurement
        [rho,~,~,rhoDot,~,~] = groundsStationOutput(Stations.(stationNames{stnID}),X(ii,:)', phi_G(ii));

        % residuals
        delta_y(nn,3) = measData(rawDataStartIndx+nn,3) - rho;
        delta_y(nn,4) = measData(rawDataStartIndx+nn,4) - rhoDot;

        % measurement sensitivity wrt current (r,v)
        Htld = calculateH_rho_rhoDot_wrt_r_v(X(ii,:)', Stations.(stationNames{stnID}), phi_G(ii));

        % STM at this time
        Phi = reshape(X(ii,7:end),6,6);
        H = Htld * Phi;    % sensitivity wrt initial state
        
        % accumulate normal eqns
        InfoMat = InfoMat + H'*(R\H);
        InfoVec = InfoVec + H'*(R\delta_y(nn,3:4)');

        % Postfit Residual
        epsilon_hat(nn,:) = delta_y(nn,3:4) - (H*delta_x_hat)';

        nn = nn + 1;
    end

    IterationData{kk}.Prefit  = delta_y(:,3:4);
    IterationData{kk}.Postfit = epsilon_hat;
    IterationData{kk}.InfoMat = InfoMat;
    IterationData{kk}.InfoVec = InfoVec;

    delta_x_hat = InfoMat \ InfoVec;   % Gauss-Newton update
    MalhonbisDist =  delta_x_hat'*(P0\delta_x_hat);
    delta_x_bar = delta_x_bar - delta_x_hat; % Update a priori state deviation

    if kk>=maxIterLimit
         warning(['Batch filter did not converge within %d iterations.\n' ...
             'Last step Mahalanobis norm: %.3e\n'], ...
             maxIterLimit, delta_x_hat'*(P0\delta_x_hat));
        break
    end
end
kk
%% Final post-fit difference from OE-based initial state
delta_x_postfit = X0(1:6) - [r0_ECI; v0_ECI];

figure(1)
figure(2)
for ii = 1:kk
    for jj = 1:2
        figure(jj)
        subplot(kk, 2, 2*ii-1)
        plot(t(indxStart:indxEnd), IterationData{ii}.Prefit(:,jj), 'o')
        title(['Prefit Residual Range', num2str(ii)])
        grid on
    
        subplot(kk, 2, 2*ii)
        plot(t(indxStart:indxEnd), IterationData{ii}.Postfit(:,jj), 'o')
        title(['Postfit Residual Range', num2str(ii)])
        grid on
    end
end
figure(1)
sgtitle('Range Pre and Postfit Residual')
figure(2)
sgtitle('Range Rate Pre and Postfit Residual')