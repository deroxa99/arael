function g = gravAcc(t,r_ijk,mu,Re,Pnm_norm,Cnm_mod,Snm_mod,n,ref_sys,obs,et)
% -------------------------------------------------------------------------
% GRAVACC - Compute Gravitational Acceleration in body-fixed Axes
% -------------------------------------------------------------------------
% Compute Gravitational Acceleration modeled by a spherical harmonic 
% expansion in body-fixed reference frame with normalized coefficients 
% (Cnm_norm,Snm_norm)
% -------------------------------------------------------------------------
% INPUTS:
% t          : [1,1] - Integration time [s]
% r          : [3,1] - position of the S/C in body-centered inertial frame
%                      [km]
% mu         : [1,1] - gravitational constant [km^3/s^2]
% Re         : [1,1] - reference radius [km]
% Pnm_norm   : [(n+1)*(n+2)/2,1] - vector of normalised associated legendre
%                                  polynomial up to order n+1
% Cnm_mod    : [(n+2)*(n+3)/2,3] - Cnm modified gravity coefficients along 
%                                  x,y,z [-]
% Snm_mod    : [(n+2)*(n+3)/2,3] - Snm modified gravity coefficients along 
%                                  x,y,z [-]
% n          : [1,1] - degree of the expansion
%  ref_sys [-]    -  String that defines the reference system as one of
%                    the following:
%                     - Earth: 'J2000'
%                     - Moon: 'MOON_PA_INERTIAL', 'MOON_ME_INERTIAL'
%
%  obs     [-]    -  String containing the name of the central body.
%                    Supported ones are:
%                     - 'EARTH'
%                     - 'MOON'
% et         : [1,1] - ephemeris time [s]
% -------------------------------------------------------------------------
% OUTPUTS:
% g          : [3,1] - Gravitational Acceleration in inertial reference 
%                      frame [km^2/s^2]
% -------------------------------------------------------------------------
% SOURCES: Construction of spherical harmonic series for the potential 
%          derivatives of arbitrary orders in the geocentric
%          Earth-fixed reference frame - M. S. Petrovskaya · A. N. Vershkov
% -------------------------------------------------------------------------
% CONTRIBUTORS: Alessio Derobertis
% -------------------------------------------------------------------------
% CHANGELOG:   2023/09/26 - V1: First draft
%              2023/10/16 - V2: Optimized by removing computation of 
%                               modified coefficients
%              2023/10/17 - V3: OPtimized by using sum functions
% -------------------------------------------------------------------------

% convert into body-fixed frame
switch obs
    case 'EARTH'
        % compute rotation matrix from ECI to ECEF
        ROT = cspice_pxform(ref_sys, 'ITRF93', t+et);

    case 'MOON'
        % compute rotation matrix from MCI to Moon Principal Axes
        ROT = cspice_pxform(ref_sys, 'MOON_PA', t+et);

    case 'MARS'
        % compute rotation matrix from MCI to Moon Principal Axes
        ROT = cspice_pxform(ref_sys, 'IAU_MARS', t+et);

    case 'VENUS'
        % compute rotation matrix from MCI to Moon Principal Axes
        ROT = cspice_pxform(ref_sys, 'IAU_VENUS', t+et);

end

r_body = ROT*r_ijk;

% Convert from cartesian to spherical
r = norm(r_body);
long = atan2(r_body(2), r_body(1));
lat = asin(r_body(3) / r);

% Evaluate the functions
Pnm_norm_eval = Pnm_norm(sin(lat));

% Initialize quantities
g = int32((n + 2) * (n + 3) / 2);
clong = zeros(n+2,1);
slong = zeros(n+2,1);
ri = zeros(n+2,1);
clong_vec = zeros(g, 1);
slong_vec = zeros(g, 1);
ri_vec = zeros(g, 1);
clong(1) = 1;
slong(1) = 0;
ri(1) = Re / r;
cl = cos(long);
sl = sin(long);

% Compute quantities
for i = 1:n + 1
    clong(i+1) = cl * clong(i) - sl * slong(i);
    slong(i+1) = cl * slong(i) + sl * clong(i);
    ri(i+1) = ri(i) * ri(1);
end

% Construct vectors
c = 1;
for i = 0:n + 1
    for j = 1:i+1
        clong_vec(c) = clong(j);
        slong_vec(c) = slong(j);
        ri_vec(c) = ri(i+1);
        c = c + 1;
    end
end

% Compute sum
sum_i = ri_vec .* Pnm_norm_eval .* (Cnm_mod .* clong_vec + Snm_mod .* slong_vec);

% Compute gravity acceleration
GMr2 = mu / Re^2;
a_x = -GMr2 * sum(sum_i(:, 1));
a_y = -GMr2 * sum(sum_i(:, 2));
a_z = GMr2 * sum(sum_i(:, 3));

% Assemble gravity vector
g_body = [a_x; a_y; a_z];

% Transform into inertial frame
g = ROT'*g_body;

end