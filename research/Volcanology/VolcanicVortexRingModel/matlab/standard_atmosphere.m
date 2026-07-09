function [T,P] = standard_atmosphere(z_asl)
% STANDARD_ATMOSPHERE  ICAO / US Standard Atmosphere 1976, T(z) and P(z).
%
%   [T,P] = standard_atmosphere(z_asl)
%
%   z_asl : geometric altitude above sea level [m], scalar or array.
%           Valid for 0-84852 m.
%   T     : temperature [K]
%   P     : pressure [Pa]
%
% Johan Gilchrist

g0 = 9.80665;      % [m/s^2]
Rd = 287.053;      % [J/(kg K)]

zb = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852];
Lb = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002, 0];

Tb = zeros(size(zb)); Pb = zeros(size(zb));
Tb(1) = 288.15; Pb(1) = 101325;
for i = 2:numel(zb)
    dz = zb(i) - zb(i-1);
    T0 = Tb(i-1); L0 = Lb(i-1); P0 = Pb(i-1);
    T1 = T0 + L0*dz;
    if abs(L0) > 1e-12
        P1 = P0*(T0/T1)^(g0/(Rd*L0));
    else
        P1 = P0*exp(-g0*dz/(Rd*T0));
    end
    Tb(i) = T1; Pb(i) = P1;
end

z_asl = z_asl(:);
T = zeros(size(z_asl));
P = zeros(size(z_asl));
for k = 1:numel(z_asl)
    i = find(zb <= z_asl(k), 1, 'last');
    if isempty(i); i = 1; end
    if i == numel(zb); i = numel(zb)-1; end
    dz = z_asl(k) - zb(i);
    T(k) = Tb(i) + Lb(i)*dz;
    if abs(Lb(i)) > 1e-12
        P(k) = Pb(i)*(Tb(i)/T(k))^(g0/(Rd*Lb(i)));
    else
        P(k) = Pb(i)*exp(-g0*dz/(Rd*Tb(i)));
    end
end

end
