function [rho,P] = air_density(T, z_asl, mode, profile, T_offset)
% AIR_DENSITY  Ideal-gas density [kg/m^3] of a dry-air-like parcel at
% temperature T [K], assumed to be in pressure equilibrium with the local
% ambient atmosphere at altitude z_asl [m a.s.l.]:
%
%   rho = P_ambient(z_asl) / (Rd * T)
%
% Generalizes the "AirDensityFunc(T,z)" convention used elsewhere in this
% project to take the ambient pressure from a full T(z),P(z) atmospheric
% profile (standard or user-supplied) instead of a fixed value.
%
%   [rho,P] = air_density(T, z_asl, mode, profile, T_offset)
%
% Johan Gilchrist

if nargin < 3; mode = 'standard'; end
if nargin < 4; profile = []; end
if nargin < 5; T_offset = 0; end

Rd = 287.053; % [J/(kg K)]
[~,P] = atmosphere_TP(z_asl, mode, profile, T_offset);
rho = P./(Rd.*T);

end
