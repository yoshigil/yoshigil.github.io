function [T,P] = atmosphere_TP(z_asl, mode, profile, T_offset)
% ATMOSPHERE_TP  Ambient temperature and pressure at altitude(s) z_asl [m a.s.l.].
%
%   [T,P] = atmosphere_TP(z_asl, mode, profile, T_offset)
%
%   z_asl    : altitude above sea level [m], scalar or array
%   mode     : 'standard' (built-in ICAO/US Standard Atmosphere 1976, default)
%              or 'profile' (interpolate a user-supplied sounding)
%   profile  : struct with fields z [m], T [K], P [Pa] (see load_atm_profile.m),
%              required if mode = 'profile'. Ignored/[] otherwise.
%   T_offset : optional uniform additive temperature offset [K] applied to
%              the standard atmosphere only (default 0) -- handy for shifting
%              the standard profile to match a known local/near-vent surface
%              temperature without needing a full custom sounding.
%
%   T [K], P [Pa]
%
% Johan Gilchrist

if nargin < 2 || isempty(mode); mode = 'standard'; end
if nargin < 3; profile = []; end
if nargin < 4 || isempty(T_offset); T_offset = 0; end

if strcmpi(mode,'profile')
    if isempty(profile) || ~isfield(profile,'z')
        error('atmosphere_TP:noProfile', ...
            ['mode=''profile'' requires a profile struct with fields z,T,P ' ...
             '(see load_atm_profile.m). Use mode=''standard'' otherwise.']);
    end
    T = interp1(profile.z, profile.T, z_asl, 'linear','extrap');
    P = interp1(profile.z, profile.P, z_asl, 'linear','extrap');
else
    [T,P] = standard_atmosphere(z_asl);
    T = T + T_offset;
end

end
