function y_modeled = G_groundStationRangeAndRangeRate(X, Stations, t, phi_G0)

phi_G = getConstants().Earth.omega*t+phi_G0;
fields = fieldnames(Stations);

rhoVec     = NaN(length(t), 1);
rhoDotVec  = NaN(length(t), 1);
IdVec      = NaN(length(t), 1);

for ii = 1:length(t)
    for jj = 1:numel(fields)
        stationName = fields{jj};
        [rho, ~, El, rhoDot, ~, ~] = groundsStationOutput(Stations.(stationName), X(ii,:)', phi_G(ii));

        if El >= 10 * getConstants().Conversions.deg2rad
            rhoVec(ii)    = rho;
            rhoDotVec(ii) = rhoDot;
            IdVec(ii)     = jj;   
            break;
        end
    end
end

y_modeled = [t(:), IdVec, rhoVec*getConstants().Conversions.m2km, rhoDotVec*getConstants().Conversions.m2km];
y_modeled = y_modeled(~any(isnan(y_modeled), 2), :);