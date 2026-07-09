function N2 = brunt_vaisala_N2(z_asl, mode, profile, T_offset, dz)
% BRUNT_VAISALA_N2  Real, height-dependent Brunt-Vaisala frequency squared
% of the ambient atmosphere, N^2(z) = (g/theta_a) d(theta_a)/dz, evaluated
% at altitude z_asl [m a.s.l.] using a centered finite difference of step
% dz [m] (default 1 m).
%
% Using potential temperature (rather than a raw density difference)
% correctly generalizes the constant-N^2 Boussinesq stratification used in
% the tank-experiment model to a real, compressible atmosphere: it
% automatically accounts for the adiabatic expansion/cooling that a rising
% ring undergoes purely from moving to lower ambient pressure, separate
% from any mixing with entrained air.
%
%   N2 = brunt_vaisala_N2(z_asl, mode, profile, T_offset, dz)
%
% Johan Gilchrist

if nargin < 2 || isempty(mode); mode = 'standard'; end
if nargin < 3; profile = []; end
if nargin < 4 || isempty(T_offset); T_offset = 0; end
if nargin < 5 || isempty(dz); dz = 1; end

g0 = 9.80665;

[Tm,Pm] = atmosphere_TP(z_asl, mode, profile, T_offset);
[T1,P1] = atmosphere_TP(z_asl-dz, mode, profile, T_offset);
[T2,P2] = atmosphere_TP(z_asl+dz, mode, profile, T_offset);

theta_m = potential_temperature(Tm,Pm);
theta1  = potential_temperature(T1,P1);
theta2  = potential_temperature(T2,P2);

dthetadz = (theta2 - theta1)./(2*dz);
N2 = g0./theta_m .* dthetadz;

end
