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
r0_ECI = r0_ECI + [0.541708589833822;0.215276743995227;-0.0246521235406893]*getConstants().Conversions.km2m;
v0_ECI = v0_ECI+[-5.6364922705741e-5;-0.000394199357037521;0.00051130212438399]*getConstants().Conversions.km2m;
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
DCO = tf;
% force ode45 output at measurement times
tspan = [t0, measData((measData(:,1)<tf & measData(:,1)>t0),1)', tf];

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
phi_G0 = 170*getConstants().Conversions.deg2rad;

%% Batch filter (Gauss–Newton on initial state)
delta_x_hat = zeros(6,1);           
delta_x_bar = zeros(6,1);
firstIteration = true;
kk = 0;
maxIterLimit = 8;
MalhonbisDist =  delta_x_hat'*(P0\delta_x_hat);
distCount = 0;
while distCount<2
    kk = kk + 1;
    X0(1:6) = X0(1:6) + delta_x_hat;      % update initial state guess
    IterationData{kk}.X0 = X0(1:6);
    IterationData{kk}.DeltaX = X0(1:6)-[r0_ECI; v0_ECI];
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
    if  MalhonbisDist<=1e-10
        distCount = distCount+1;
    else
        distCount = 0;
    end
    delta_x_bar = delta_x_bar - delta_x_hat; % Update a priori state deviation

    if kk>=maxIterLimit
         warning(['Batch filter did not converge within %d iterations.\n' ...
             'Last step Mahalanobis norm: %.3e\n'], ...
             maxIterLimit, delta_x_hat'*(P0\delta_x_hat));
        break
    end
end
kk
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

%% 36 hour truth Data
% Initial OE 
OE0.a     = 10000*getConstants().Conversions.km2m;
OE0.e     = 0.001;
OE0.i     = 40*getConstants().Conversions.deg2rad;
OE0.Omega = 80*getConstants().Conversions.deg2rad;
OE0.omega = 40*getConstants().Conversions.deg2rad;
OE0.theta = 0;
[r0_ECI, v0_ECI] = PositionVelocityFromOrbitalElements(OE0);

% Build augmented initial state (r,v,STM)
X0 = [r0_ECI; v0_ECI; reshape(eye(6),36,1)];
t0 = 0;
tf = 36*getConstants().Conversions.hr2s;

% force ode45 output at measurement times
tspan = t0:10:tf;

phi_G0 = 170*getConstants().Conversions.deg2rad;
[t,X_truth] = ode45(@twobodyEOMJ2WithSTM, tspan, X0, options);
[t,X_est] = ode45(@twobodyEOMJ2WithSTM, tspan, [IterationData{end}.X0;reshape(eye(6),36,1)], options);

sigma_r = nan(length(t), 3);   % [σ_rx, σ_ry, σ_rz]
sigma_v = nan(length(t), 3);   % [σ_vx, σ_vy, σ_vz]

for ii = 1:length(t)
    Phi = reshape(X_est(ii,7:end), 6, 6);      % STM at time t(k)
    Pk  = Phi * (IterationData{end}.InfoMat\Phi');            % covariance at time t(k)

    % 1σ position (x,y,z)
    sigma_r(ii,:) = sqrt(diag(Pk(1:3,1:3))).';

    % 1σ velocity (vx,vy,vz)
    sigma_v(ii,:) = sqrt(diag(Pk(4:6,4:6))).';

    R = X_est(ii,1:3)/norm(X_est(ii,1:3));
    N = cross(X_est(ii,1:3),X_est(ii,4:6))/norm(cross(X_est(ii,1:3),X_est(ii,4:6)));
    T = cross(N,R);
    DCM_RTN2ECI = [R,T,N];
    DCM_ECI2RTN = DCM_RTN2ECI';
    X_est_RTN(ii,1:3) = DCM_ECI2RTN*X_est(ii,1:3);
    X_est_RTN(ii,4:6) = DCM_ECI2RTN*X_est(ii,4:6);
    X_truth_RTN(ii,1:3) = DCM_ECI2RTN*X_truth(ii,1:3);
    X_truth_RTN(ii,4:6) = DCM_ECI2RTN*X_truth(ii,4:6);
    Pk_RTN = [DCM_ECI2RTN, zeros(3);zeros(3),DCM_ECI2RTN] * Pk *  [DCM_ECI2RTN, zeros(3);zeros(3),DCM_ECI2RTN]';
    % 1σ position (x,y,z)
    sigma_r_RTN(ii,:) = sqrt(diag(Pk_RTN(1:3,1:3))).';

    % 1σ velocity (vx,vy,vz)
    sigma_v_RTN(ii,:) = sqrt(diag(Pk_RTN(4:6,4:6))).';
end

figure
ax(1) = subplot(3,1,1);
plot(t*getConstants().Conversions.s2hrs, X_est(:,1) - X_truth(:,1))
hold on
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
plot(t*getConstants().Conversions.s2hrs, sigma_r(:,1)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_r(:,1)*3, '--k')
legend('State_{diff}', 'DC0', '+/-3 \sigma bound')
ylabel('r_{x,ECI} [m]')
grid on
ax(2) = subplot(3,1,2);
plot(t*getConstants().Conversions.s2hrs, X_est(:,2) - X_truth(:,2))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_r(:,2)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_r(:,2)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('r_{y,ECI} [m]')
grid on
ax(3) = subplot(3,1,3);
plot(t*getConstants().Conversions.s2hrs, X_est(:,3) - X_truth(:,3))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_r(:,3)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_r(:,3)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('r_{z,ECI} [m]')
grid on
xlabel('Time [Hours]')
linkaxes(ax, 'x');
sgtitle('Position Difference in Estimate and Truth')

figure
ax(1) = subplot(3,1,1);
plot(t*getConstants().Conversions.s2hrs, X_est(:,4) - X_truth(:,4))
hold on
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
plot(t*getConstants().Conversions.s2hrs, sigma_v(:,1)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_v(:,1)*3, '--k')
ylabel('v_{x,ECI} [m/s]')
grid on
legend('State_{diff}', 'DC0', '+/-3 \sigma bound')
ax(2) = subplot(3,1,2);
plot(t*getConstants().Conversions.s2hrs, X_est(:,5) - X_truth(:,5))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_v(:,2)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_v(:,2)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('v_{y,ECI} [m/s]')
grid on
ax(3) = subplot(3,1,3);
plot(t*getConstants().Conversions.s2hrs, X_est(:,6) - X_truth(:,6))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_v(:,3)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_v(:,3)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('v_{z,ECI} [m/s]')
grid on
xlabel('Time [Hours]')
linkaxes(ax, 'x');
sgtitle('Velocity Difference in Estimate and Truth')

figure
ax(1) = subplot(3,1,1);
plot(t*getConstants().Conversions.s2hrs, X_est_RTN(:,1) - X_truth_RTN(:,1))
hold on
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
plot(t*getConstants().Conversions.s2hrs, sigma_r_RTN(:,1)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_r_RTN(:,1)*3, '--k')
legend('State_{diff}', 'DC0', '+/-3 \sigma bound')
ylabel('r_{x,ECI} [m]')
grid on
ax(2) = subplot(3,1,2);
plot(t*getConstants().Conversions.s2hrs, X_est_RTN(:,2) - X_truth_RTN(:,2))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_r_RTN(:,2)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_r_RTN(:,2)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('r_{y,ECI} [m]')
grid on
ax(3) = subplot(3,1,3);
plot(t*getConstants().Conversions.s2hrs, X_est_RTN(:,3) - X_truth_RTN(:,3))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_r_RTN(:,3)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_r_RTN(:,3)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('r_{z,ECI} [m]')
grid on
xlabel('Time [Hours]')
linkaxes(ax, 'x');
sgtitle('Position Difference in Estimate and Truth')

figure
ax(1) = subplot(3,1,1);
plot(t*getConstants().Conversions.s2hrs, X_est_RTN(:,4) - X_truth_RTN(:,4))
hold on
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
plot(t*getConstants().Conversions.s2hrs, sigma_v_RTN(:,1)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_v_RTN(:,1)*3, '--k')
ylabel('v_{x,ECI} [m/s]')
grid on
legend('State_{diff}', 'DC0', '+/-3 \sigma bound')
ax(2) = subplot(3,1,2);
plot(t*getConstants().Conversions.s2hrs, X_est_RTN(:,5) - X_truth_RTN(:,5))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_v_RTN(:,2)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_v_RTN(:,2)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('v_{y,ECI} [m/s]')
grid on
ax(3) = subplot(3,1,3);
plot(t*getConstants().Conversions.s2hrs, X_est_RTN(:,6) - X_truth_RTN(:,6))
hold on
plot(t*getConstants().Conversions.s2hrs, sigma_v_RTN(:,3)*3, '--k')
plot(t*getConstants().Conversions.s2hrs, -sigma_v_RTN(:,3)*3, '--k')
xline(DCO*getConstants().Conversions.s2hrs, 'r--')
ylabel('v_{z,ECI} [m/s]')
grid on
xlabel('Time [Hours]')
linkaxes(ax, 'x');
sgtitle('Velocity Difference in Estimate and Truth')