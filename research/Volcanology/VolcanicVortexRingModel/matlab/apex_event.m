function [value,isterminal,direction] = apex_event(~, uvF, u_min)
% APEX_EVENT  ode23s Events function that stops integration when the
% ring's vertical velocity u = M/V decays to u_min [m/s].
%
% The height-domain momentum equation dM/dz = (2/3) B V / M is singular as
% u -> 0 (the ring's rise apex), so integration should be stopped at/near
% that point rather than letting the solver grind through ever-shrinking
% step sizes trying to push past it. The physically meaningful "rise
% height" is the level of neutral buoyancy (where g' crosses zero), which
% is reached, and well resolved, before this singularity.
%
% Johan Gilchrist

if nargin < 3 || isempty(u_min); u_min = 1e-3; end
V = uvF(1); M = uvF(2);
value = M/V - u_min;
isterminal = 1;
direction = -1;

end
