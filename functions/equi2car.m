function rv = equi2car(equi, mu)
% -------------------------------------------------------------------------
% EQUI2CAR - Convert equinoctial elements to cartesian ones
% -------------------------------------------------------------------------
% Convert a set of equinoctial elements into a cartesian vector containing 
% position and velocity expressed into the intertial basis
% -------------------------------------------------------------------------
% INPUTS:
% equi         : [6,1] - equinoctial elements of the orbiter. Ordered in 
%                        the following way: [p, f, g, h, k, L]
% mu           : [1,1] - gravitational constant of the attractor [km^3/s^2]
% -------------------------------------------------------------------------
% OUTPUTS:
% rv           : [6,1] - cartesian state vector in inertial reference
%                        frame [km,km/s]
% -------------------------------------------------------------------------
% CONTRIBUTORS: Alessio Derobertis
% -------------------------------------------------------------------------

% unpack equinoctial elements
p = equi(1);
f = equi(2);
g = equi(3);
h = equi(4);
k = equi(5);
L = equi(6);

% abbreviations
sL = sin(L);
cL = cos(L);
s2 = 1 + h^2 + k^2;
a2 = h^2 - k^2;
w = 1 + f * cL + g * sL;
r = p / w;

% initialize state vector
rv = zeros(6, 1);
rv(1:3) = (r / s2) * [cL + a2 * cL + 2 * h * k * sL;
    sL - a2 * sL + 2 * h * k * cL;
    2 * (h * sL - k * cL)];
rv(4:6) = sqrt(mu / p) * (1 / s2) * [-(sL + a2 * sL - 2 * h * k * cL + g - 2 * f * h * k + a2 * g);
    -(-cL + a2 * cL + 2 * h * k * sL - f + 2 * g * h * k + a2 * f);
    2 * (h * cL + k * sL + f * h + g * k)];
end