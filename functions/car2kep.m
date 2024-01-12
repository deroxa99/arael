function kep = car2kep(rv, mu)
% ------------------------------------------------------------------------
% CAR2KEP - Convert cartesian inertial state into keplerian elements
% ------------------------------------------------------------------------
% Convert a cartesian vector containing position and velocity expressed in 
% the inertial basis into a set of keplerian elements
% ------------------------------------------------------------------------
% INPUTS:
% rv          : [6,1] - cartesian state vector in inertial reference frame
% mu          : [1,1] - gravitational parameter of the attractor [km^3/s^2]
%
% OUTPUTS:
% kep         : [6,1] - keplerian elements of the orbiter. Ordered in the 
%                       following way:
%                       - SMA:    SEMI-MAJOR AXIS [km]
%                       - ECC:    ECCENTRICITY [-]
%                       - INC:    INCLINATION [rad]
%                       - LNODE:  LONGITUDE OF THE ASCENDING NODE [rad]
%                       - ARGP:   ARGUMENT OF PERIAPSIS [rad]
%                       - THETA:  TRUE ANOMALY [rad]
% ------------------------------------------------------------------------
% CONTRIBUTORS: Alessio Derobertis
% ------------------------------------------------------------------------

% unpack vector
r = rv(1:3);
v = rv(4:6);

% Position Vector and Velocity Vector norms
rmod = norm(r);
vmod = norm(v);

% Specific Angular Momentum
hvec = cross(r, v);
h = norm(hvec);

% Inclination
INC = acos(hvec(3)/h);

% Eccentricity Vector and Eccentricity
evec = (1/mu) * ((vmod^2 - mu/rmod)*r - dot(r,v)*v);
ECC = norm(evec);

% Specific Mechanical Energy and Semi-major Axis
eps = 0.5*(vmod^2) - mu/rmod;
SMA = - mu/(2*eps);

% Axis of Nodes
K = [0; 0; 1];
Nvec = cross(K, hvec);
N = norm(Nvec);

% Right Ascension of the Ascending Node
if Nvec(2) >= 0
    LNODE = acos(Nvec(1)/N);
else
    LNODE = 2*pi - acos(Nvec(1)/N);
end

% Pericentre Anomaly
if evec(3) >= 0
    ARGP = acos(dot(Nvec, evec)/(N*ECC));
else
    ARGP = 2*pi - acos(dot(Nvec, evec)/(N*ECC));
end

% Radial Velocity
vr = dot(r, v)/rmod;

% True Anomaly
if vr >= 0
    THETA = acos(dot(evec, r)/(ECC*rmod));
else
    THETA = 2*pi - acos(dot(evec, r)/(ECC*rmod));
end

% assemble keplerians
kep = [SMA; ECC; INC; LNODE; ARGP; THETA];

end