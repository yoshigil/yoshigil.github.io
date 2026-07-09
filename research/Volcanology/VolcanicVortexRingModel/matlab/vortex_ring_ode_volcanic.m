function dz = vortex_ring_ode_volcanic(z_above_vent, uvF, p)
% VORTEX_RING_ODE_VOLCANIC  Height-domain ODE right-hand side for the 1D
% variable-entrainment buoyant vortex ring / thermal rise model,
% generalized to a real (compressible) atmosphere.
%
% Based on models from Morton et al. 1956; Turner 1962; Sparks et al. 1997.
% This is the volcanic-atmosphere counterpart of
% thermalfuncup_height_var_alpha.m.
%
% State vector uvF = [V; M; B; alpha]
%   V     : ring volume                         [m^3]      ~ (4/3) pi r^3
%   M     : ring volume x vertical velocity      [m^4/s]    ~ V*u
%   B     : ring volume x reduced gravity        [m^4/s^2]  ~ V*g'
%   alpha : entrainment coefficient (diagnostic state, tracked so it can be
%           output alongside V,M,B)
%
%   dV/dz = 3 * alpha * V^(2/3)
%   dM/dz = (2/3) * B * V / M
%   dB/dz = -N^2(z) * V
%
% The only physics that changes relative to the analog tank-experiment
% model is:
%   (1) B is defined via potential-temperature excess rather than a
%       fixed-fluid density excess, so the ring's own adiabatic
%       expansion/cooling as it rises to lower ambient pressure is
%       captured automatically; and
%   (2) N^2(z) is the REAL, height-varying stratification of the
%       atmosphere (standard atmosphere or a supplied sounding),
%       evaluated at the true altitude p.z0_asl + z_above_vent.
%
% p : parameter struct with fields
%       z0_asl     - vent altitude above sea level [m]
%       atmMode    - 'standard' or 'profile'
%       atmProfile - struct (or []) for atmosphere_TP.m
%       T_offset   - [K] standard-atmosphere temperature offset
%       c          - entrainment closure coefficient (Turner 1962)
%       alpha_max  - hard cap on the entrainment coefficient (see below)
%
% alpha_max: a volcanic ring exiting into cold, high-altitude ambient air
% can have a MUCH larger relative buoyancy anomaly than the small-
% density-contrast analog tank experiments the alpha(Ri) closure was built
% for (fractional potential-temperature excess of order 1 rather than a
% few percent). Left uncapped, the closure's positive feedback (higher
% alpha -> more dilution of momentum -> even higher alpha) can produce a
% genuine finite-height blow-up of V. Capping alpha is both a standard
% numerical safeguard for variable-entrainment closures of this type (see
% e.g. Kaminski et al. 2005) and physically grounded here: fitted near-vent
% ("jet development") entrainment coefficients for Etna's Rings vent
% events reach up to ~0.95 (see INGV/Etna_Rings_Vent/Rings_alphae.m), so
% alpha_max = 1.0 is a reasonable, data-informed ceiling.
%
% Johan Gilchrist

dz = zeros(4,1);
V = uvF(1); M = uvF(2); B = uvF(3);

r = (3*V/(4*pi))^(1/3);
u = M/V;

alpha = abs(B)/(2*pi*p.c*(r*u)^2);
alpha = min(alpha, p.alpha_max);

N2 = brunt_vaisala_N2(p.z0_asl + z_above_vent, p.atmMode, p.atmProfile, p.T_offset);

dVdz = 3*alpha*V^(2/3);
dMdz = (2/3)*B*V/M;
dBdz = -N2*V;

dz(1) = dVdz;
dz(2) = dMdz;
dz(3) = dBdz;
dz(4) = alpha - uvF(4);

end
