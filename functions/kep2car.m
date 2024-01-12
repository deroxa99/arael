function [x] = kep2car(kep,mu)
% ------------------------------------------------------------------------
% KEP2CAR - Conversion from Keplerian coordinates to Cartesian elements. 
% Angles in radians.
% ------------------------------------------------------------------------
%
% INPUT:
%   kep           [6X1]   Vector containing the keplerian elements of the
%                         body, in this order: 
%                         a - Semi-major Axis                         [km]
%                         e - Eccentricity                            [-]
%                         i - Inclination                             [rad]
%                         RAAN - Right Ascesion of Ascending Node     [rad]
%                         omega - Pericentre Anomaly                  [rad]
%                         f0 - True Anomaly                           [rad]
%   mu            [1x1]   Gravitational parameter                [km^3/s^2]
% 
% OUTPUT:
%   r             [3x1]   Position vector                    [km]
%   v             [3x1]   Velocity vector                    [km/s]
% ------------------------------------------------------------------------
% Contributors:
%   Matteo Gant, Pasquale Castellano, Giovanni Fereoli, Alessio Derobertis
% ------------------------------------------------------------------------
% Version:
%   V2: 2022-23-03 || Updated - Alessio Deroberts
%                   
%   V1: 2021-12-31 || Matteo Gant, Pasquale Castellano, Giovanni Fereoli, 
%                     Alessio Derobertis
% ------------------------------------------------------------------------

% unpack the keplerian
a = kep(1);
e = kep(2);
i = kep(3);
RAAN = kep(4);
omega = kep(5);
f0 = kep(6);

p = a*(1 - e^2);  % Semi-Latus Rectum [km]
r_norm = p/(1 + e*cos(f0));  % Position Modulus [km]

% State Vectors in Perifocal Frame
r_PF = r_norm*[cos(f0); sin(f0); 0];        % Position Vector [km]
v_PF = sqrt(mu/p)*[-sin(f0); e + cos(f0); 0]; % Velocity vector [km/s] 

% Rotation Matrix around axis 3 of angle OMEGA  = 
R3_OMEGA = [cos(RAAN), sin(RAAN), 0; -sin(RAAN), cos(RAAN) , 0; 0, 0, 1];
% Rotation Matrix around axis 1 of angle i = Inclination 
R1_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
% Rotation Matrix around axis 3 of angle omega = 
R3_omega = [cos(omega), sin(omega), 0; -sin(omega), cos(omega) ,0 ; 0, 0, 1];

% Overall Rotation Matrix
% PF Perifocal Frame --> ECEI Earth-Centred Equatorial Inertial
R_ECEI2PF = R3_omega*R1_i*R3_OMEGA;
R_PF2ECEI = R_ECEI2PF';

% State Vectors in ECEI Earth-Centred Equatorial Inertial
r = R_PF2ECEI*r_PF;  % Position Vector [km]
v = R_PF2ECEI*v_PF;  % Velocity vector [km/s] 

x = [r;v];
end