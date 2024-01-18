function [aSRP] = pertSRP_full(t,rSC,set,cR,A,m,ref_sys,obs,et)
% ------------------------------------------------------------------------
% DESCRIPTION:
% AVG3B - Compute averaged perturbing acceleration due to solar radiation
%         pressure. Specialized for 'full' mode.
% ------------------------------------------------------------------------
% INPUT ARGUMENTS:
%  t       [1,1]  -  Integration time [s]
%
%  rSC     [3,1]  -  Position of the S/C in inertial frame [km]
%
%  set     [-]    - Set as:
%                   - 'on': SRP is included
%                   - 'off': SRP is not included
%
%  cR      [1,1]  -  Reflectivity coefficient of the spacecraft [-]
%
%  A       [1,1]  -  Surface of the spacecraft [m^2]
%
%  m       [1,1]  -  Mass of the spacecraft [kg]
%
%  ref_sys [-]    -  String that defines the reference system as one of
%                    the following:
%                     - Earth: 'J2000'
%                     - Moon: 'MOON_PA_INERTIAL', 'MOON_ME_INERTIAL'
%
%  obs     [-]    -  String containing the name of the central body.
%                    Supported ones are:
%                     - 'EARTH' (n_max = 6)
%                     - 'MOON'  (n_max = 2)
%
%  et      [1,1]  -  Initial time in seconds past J2000 
%
% OUTPUT ARGUMENTS:
%  a3B        [3x1]:   Perturbing acceleration in body-centered interial 
%                      frame [km/s]
% ------------------------------------------------------------------------
% CONTRIBUTOS: 
%  Alessio Derobertis
% ------------------------------------------------------------------------
% REFERENCE: Fundamentals of Astrodynamics and Applications - David A.
%            Vallado - Fourth edition - pag. 299/578
% ------------------------------------------------------------------------
% CHANGELOG:
%  18/01/2024 - First draft - ALessio Derobertis
% ------------------------------------------------------------------------

switch set
    case 'on'

        % constants
        AU = 149597870.7; % astronomic unit [km]
        p = 4.57e-06; % average solar radiation pressure @ 1AU [N/m^2]

        % compute position of the sun
        rSun = cspice_spkpos('SUN',t+et,ref_sys,'NONE',obs);
        rRel = rSun - rSC;
        rRelNorm = norm(rRel);

        % compute solar radiation pressure
        aSRP = -(p*cR*A/m)*(rRel/rRelNorm)*(AU^2/rRelNorm^2);

    case 'off'

        % no srp
        aSRP = [0;0;0];

end

end