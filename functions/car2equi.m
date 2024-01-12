function equi = car2equi(rv, mu)
% -------------------------------------------------------------------------
% CAR2KEP - Convert cartestian inertial state into equinoctial element
% -------------------------------------------------------------------------
% Convert a cartesian vector containing position and velocity expressed in
% the intertial basis into a set of equinoctial element
% -------------------------------------------------------------------------
% INPUTS
% rv           : [6,1] - cartesian state vector in inertial reference frame
%                        [km,km/s]
% mu           : [1,1] - gravitational parameter of the attractor
%                        [km^3/s^2]
% -------------------------------------------------------------------------
% OUTPUTS
% kep          : [6,1] - equinoctial elements of the orbiter. Ordered in
%                        the following way: [p, f, g, h, k, L]
% -------------------------------------------------------------------------
% CONTRIBUTORS: Alessio Derobertis
% -------------------------------------------------------------------------

% unpack the state
r = rv(1:3);
v = rv(4:6);

% compute equinoctial elements
rdv = dot(r, v);
rhat = r / norm(r);
rmag = norm(r);
hvec = cross(r, v);
hmag = norm(hvec);
hhat = hvec / hmag;
vhat = (rmag * v - rdv * rhat) / hmag;
p = hmag^2 / mu;
k = hhat(1) / (1.0 + hhat(3));
h = -hhat(2) / (1.0 + hhat(3));
kk = k^2;
hh = h^2;
s2 = 1.0 + hh + kk;
tkh = 2.0 * k * h;
ecc = cross(v, hvec) / mu - rhat;

fhat = [1.0 - kk + hh, tkh, -2.0 * k] / s2;
ghat = [tkh, 1.0 + kk - hh, 2.0 * h] / s2;
f = dot(ecc, fhat);
g = dot(ecc, ghat);
L = atan2(rhat(2) - vhat(1), rhat(1) + vhat(2));

% assemble equinoctial elements
equi = [p, f, g, h, k, L];

end