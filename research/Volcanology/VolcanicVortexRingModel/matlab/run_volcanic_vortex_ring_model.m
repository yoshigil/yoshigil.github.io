%% 1D Variable-Entrainment Vortex Ring / Thermal Rise Model
%% -- VOLCANIC ATMOSPHERE VERSION --
%
% Adapted from the analog buoyant-vortex-ring (thermal) tank-experiment
% model (OneDThermalModel_PureFluid_height.m / thermalfuncup_height_var_alpha.m)
% for NATURAL VOLCANIC vortex rings rising through a real, stratified
% atmosphere.
%
% The governing entrainment/momentum/buoyancy ODEs are unchanged from the
% tank-experiment model; what changes is the physics of the ambient medium
% and the source:
%
%   - The working fluid is air (ideal gas): density depends on temperature
%     AND on the local ambient pressure (which decreases with altitude),
%     not on a fixed lab-fluid density.
%   - Buoyancy is tracked via potential temperature so that the model
%     automatically captures the adiabatic expansion/cooling the ring
%     undergoes as it rises to lower ambient pressure (see
%     brunt_vaisala_N2.m).
%   - The ambient stratification N^2(z) is computed from a real
%     atmospheric T(z), P(z) profile -- either the built-in ICAO/US
%     Standard Atmosphere 1976 (standard_atmosphere.m), or a user-supplied
%     sounding (load_atm_profile.m) -- evaluated at the TRUE altitude
%     (vent altitude + height risen).
%   - The vent (source) can be placed at any altitude above sea level.
%   - No particle/ash loading is tracked (pure-gas ring), matching the
%     "PureFluid" scope of the original model and the assumption that
%     ring density is that of air at the ring's temperature.
%
% Source condition defaults below are for Etna's "Rings" vent (Bocca Nuova
% crater, ~3380 m a.s.l.), consistent with thermal-camera measurements of
% vortex ring events described elsewhere in this project (see
% INGV/Etna_Rings_Vent/EtnaRings_diam_alpha_vel_250902.m). Exit temperature
% T0 = 450 degC is a measured value, consistent with the ~280-545 degC
% range of peak temperatures observed across 11 recorded Rings events.
% Vent radius and exit velocity are PLACEHOLDERS -- replace with the
% measured values for the specific event you want to model (observed
% ranges in the same dataset: ring radius ~0-19 m, vertical velocity
% ~0.05-33 m/s).
%
% Johan Gilchrist

clear
close all

%% ------------------- USER-DEFINED PARAMETERS -------------------------
% --- Entrainment / model parameters ---
alphaup = 0.18;    % top-hat entrainment coefficient, used only for the
                    % plume-scaling comparison law below (not in the
                    % variable-alpha ODE)
c       = 1;        % O(1) coefficient in the variable-alpha closure
                    % (matches the original tank-model convention in
                    % thermalfuncup_height_var_alpha.m;
                    % EtnaRings_diam_alpha_vel_250902.m instead used
                    % c=0.26, Turner 1962 -- change here if desired).
alpha_max = 1.0;    % hard cap on the variable entrainment coefficient --
                    % see the detailed comment in
                    % vortex_ring_ode_volcanic.m. Needed because volcanic
                    % source/ambient temperature contrasts are far larger
                    % than the tank experiments the alpha(Ri) closure was
                    % tuned on; without a cap the model can blow up to NaN
                    % within a few source radii of the vent. 1.0 is
                    % consistent with the largest near-vent entrainment
                    % coefficients fitted to real Rings-vent events (~0.95).

% --- Source (vent) conditions -- EDIT THESE FOR YOUR EVENT -------------
z0_asl = 3380;      % [m] vent altitude above sea level.
                    % Default: Etna's Bocca Nuova crater ("Rings" vent).
r0     = 5;          % [m] vent/initial ring radius -- PLACEHOLDER, edit to
                    % the measured value for the event being modeled
                    % (observed range across 11 Rings events: ~0-19 m).
u0     = 20;          % [m/s] initial (exit) vertical velocity -- measured
                    % value for the event being modeled (observed range:
                    % ~0.05-33 m/s).
Tj0_C  = 450;        % [degC] initial ring temperature at the vent --
                    % measured value for Etna's "Rings" vent.
Tj0    = Tj0_C + 273.15; % [K]

% --- Atmosphere ---------------------------------------------------------
atmMode = 'standard'; % 'standard' = built-in ICAO/US Standard Atmosphere 1976
                       % 'profile'  = user-supplied T(z), P(z) sounding
atmProfile = [];

Ta0_target_C = 15;     % [degC] measured/desired ambient temperature AT THE
                       % VENT ALTITUDE. The standard atmosphere alone is
                       % generally too cold near a volcanic summit (local
                       % heating, etc.), so T_offset below is solved for
                       % automatically to make atmosphere_TP(z0_asl,...)
                       % return exactly this value; the rest of the
                       % standard-atmosphere profile is shifted by the same
                       % constant amount. Ignored if atmMode = 'profile'
                       % (supply a real sounding instead in that case).
[Ta0_std,~] = atmosphere_TP(z0_asl, atmMode, atmProfile, 0);
T_offset = (Ta0_target_C + 273.15) - Ta0_std; % [K]

% Example: use a measured/radiosonde profile instead of the standard atmosphere
% atmProfile = load_atm_profile('etna_sounding.csv','TempUnits','C','PressUnits','hPa');
% atmMode = 'profile';

% --- Integration range ---------------------------------------------------
zmax_rise = 3000;   % [m] max height ABOVE THE VENT to integrate to
zmin_rise = 2*r0;    % [m] start a couple of vent-radii above the vent
                     % (avoids the r=0 singularity in the entrainment law)
n_pts = 10000;
z = linspace(zmin_rise, zmax_rise, n_pts);

%% ------------------- INITIAL CONDITIONS -------------------------
[Ta0,Pa0] = atmosphere_TP(z0_asl, atmMode, atmProfile, T_offset);
theta_a0  = potential_temperature(Ta0,Pa0);
theta_j0  = potential_temperature(Tj0,Pa0); % ring assumed at local ambient P

V0 = (4/3)*pi*r0^3;      % [m^3] initial ring volume
M0 = V0*u0;               % [m^4/s] initial momentum
B0 = V0*9.80665*(theta_j0-theta_a0)/theta_a0; % [m^4/s^2] initial buoyancy flux
K0 = u0*r0;                % [m^2/s] initial circulation scale
gprime0 = B0/V0; % [m/s^2] initial reduced gravity (potential-temperature based)

% Source (initial) Richardson number, using the SAME reduced-gravity
% definition as the rest of the model (potential temperature -- see
% brunt_vaisala_N2.m) and as the Ri_z height profile below, so this is
% exactly Ri_z's starting value and is consistent with alpha_k.
Ri0 = abs(gprime0)*r0/u0^2;

alpha_k = min(abs(B0)/(2*pi*c*K0^2), alpha_max);

uvF = [V0; M0; B0; alpha_k]; % state vector: [V; M; B; alpha]

%% ------------------- INTEGRATE -------------------------
% NOTE: as with the original tank-experiment model, this height-domain
% formulation (dM/dz = (2/3) B V / M) is mathematically singular as u -> 0,
% i.e. right at the top of the ring's rise / near the level of neutral
% buoyancy. The apex_event below stops integration cleanly at that point
% rather than letting the solver grind through ever-shrinking step sizes.
% The physically meaningful "rise height" is the level of neutral buoyancy
% (z_LNB below), which is reached, and well resolved, before that point.

p.z0_asl = z0_asl; p.atmMode = atmMode; p.atmProfile = atmProfile;
p.T_offset = T_offset; p.c = c; p.alpha_max = alpha_max;

opts = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
    'Events',@(zz,uv) apex_event(zz,uv,1e-3));
[z_out, dz, ze, ye, ie] = ode23s(@(zz,uv) vortex_ring_ode_volcanic(zz,uv,p), z, uvF, opts);

if ~isempty(ie)
    fprintf(['[info] ring reached its rise apex (u -> 0) at z = %.1f m ' ...
        'above the vent, before the requested zmax_rise = %.0f m.\n'], ...
        z_out(end), zmax_rise);
end

V = dz(:,1); M = dz(:,2); B = dz(:,3); alpha_k_z = dz(:,4);
u = M./V;
r = (3*V/(4*pi)).^(1/3);
gprime = B./V; % [m/s^2] reduced gravity (potential-temperature based)
Ri_z = abs(gprime).*r./u.^2; % local Richardson number height profile

% Recover ring / ambient temperature profiles for diagnostics
[Ta_z,Pa_z] = atmosphere_TP(z0_asl+z_out, atmMode, atmProfile, T_offset);
theta_a_z = potential_temperature(Ta_z,Pa_z);
theta_p_z = theta_a_z + gprime.*theta_a_z./9.80665;
Rd = 287.053; cp = 1004.7; kappa = Rd/cp; P0ref = 100000;
Tp_z = theta_p_z.*(Pa_z./P0ref).^kappa;

%% ------------------- SCALING LAWS (approximate, evaluated at vent) -------------------------
% These classical scalings assume constant N^2 and are only leading-order
% estimates in a real, height-varying atmosphere.
Cth = 1; Cpl = 2.3;
N2_vent = brunt_vaisala_N2(z0_asl, atmMode, atmProfile, T_offset);

zmax_thermal_scaling  = Cth*alpha_k^(-3/4)*B0^(1/4)*N2_vent^(-1/4) - r0/alpha_k;
zmax_vortring_scaling = Cth*K0^(3/2)*B0^(-1/2)*N2_vent^(-1/4);
zmax_plume_scaling    = Cpl*alphaup^(-1/2)*B0^(1/4)*N2_vent^(-3/8);

% Level of neutral buoyancy
[~,ilnb] = min(abs(gprime));
z_LNB = z_out(ilnb);

% Height-averaged entrainment coefficient over the full integrated rise
alpha_k_z_avg = mean(alpha_k_z);
% ... and over just the rise up to the level of neutral buoyancy
alpha_k_z_avg_to_LNB = mean(alpha_k_z(1:ilnb));

fprintf('Vent altitude:                     %8.1f m a.s.l.\n', z0_asl);
fprintf('Ambient T, P at vent:               %6.1f degC, %7.1f hPa\n', Ta0-273.15, Pa0/100);
fprintf('Initial ring T, radius, velocity:   %6.1f degC, %5.2f m, %5.2f m/s\n', Tj0_C, r0, u0);
fprintf('Source Richardson number Ri0:       %8.4f\n', Ri0);
fprintf('Initial entrainment coeff. alpha_k: %8.4f\n', alpha_k);
fprintf('N^2 at vent:                        %8.3e s^-2\n', N2_vent);
fprintf('--- classical scaling estimates (approx, const-N^2) ---\n');
fprintf('  thermal scaling:   %8.1f m\n', zmax_thermal_scaling);
fprintf('  vortex-ring scal.: %8.1f m\n', zmax_vortring_scaling);
fprintf('  plume scaling:     %8.1f m\n', zmax_plume_scaling);
fprintf('--- integrated model result ---\n');
fprintf('Max height risen above vent:        %8.1f m\n', z_out(end));
fprintf('Level of neutral buoyancy (rise):    %8.1f m\n', z_LNB);
fprintf('Level of neutral buoyancy (a.s.l.):  %8.1f m\n', z0_asl+z_LNB);
fprintf('Height-averaged alpha_e (0-LNB):    %8.4f\n', alpha_k_z_avg_to_LNB);
fprintf('Height-averaged alpha_e (full rise): %8.4f\n', alpha_k_z_avg);

%% ------------------- PLOTS -------------------------
% Figure 1: fluxes (Volume, Momentum, Buoyancy)
figure;
    subplot(131); plot(V,z_out,'.')
        title('Volume'); ylabel('Height above vent [m]'); xlabel('V [m^3]')
    subplot(132); plot(M,z_out,'.')
        title('Momentum'); xlabel('M [m^4/s]')
    subplot(133); plot(B,z_out,'.')
        title('Buoyancy'); xlabel('B [m^4/s^2]')

% Figure 2: profiles (velocity, radius, reduced gravity, temperature)
figure;
    subplot(141); plot(u,z_out,'.'); xlabel('Vertical velocity [m/s]')
        ylabel('Height above vent [m]')
    subplot(142); plot(r,z_out,'.'); xlabel('Radius [m]')
    subplot(143); plot(gprime,z_out,'.'); xlabel('g_{prime} [m/s^2]')
    subplot(144); plot(Tp_z-273.15,z_out,'.'); hold on
        plot(Ta_z-273.15,z_out,'--')
        xlabel('Temperature [deg C]'); legend('ring','ambient')

% Figure 3: entrainment coefficient and Richardson number. Ri = |g'| r /
% u^2 diverges as u -> 0 near the rise apex, so (as in the original
% tank-experiment model) it's plotted on a log x-axis.
figure;
    subplot(121); plot(alpha_k_z,z_out,'.')
        xlabel('\alpha_e'); ylabel('Height above vent [m]'); xlim([0 1])
    subplot(122); semilogx(Ri_z,z_out,'.')
        xlabel('Richardson number, Ri')
