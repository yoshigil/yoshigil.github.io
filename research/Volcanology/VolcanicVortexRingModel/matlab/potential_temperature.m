function theta = potential_temperature(T,P)
% POTENTIAL_TEMPERATURE  Dry potential temperature theta = T*(P0ref/P)^kappa.
%
% Johan Gilchrist

Rd = 287.053; cp = 1004.7; kappa = Rd/cp; P0ref = 100000; % [Pa]
theta = T.*(P0ref./P).^kappa;

end
