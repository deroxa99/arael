function [aGF_iner] = pertZH(t,r_iner,ref_sys,obs,mu,R,n,JN,et)
% ------------------------------------------------------------------------
% DESCRIPTION:
% pertZH - Compute perturbing acceleration due to gravity zonal harmonics
% ------------------------------------------------------------------------
% INPUT ARGUMENTS:
%  t       [1,1]  -  Integration time [s]
%
%  rSC     [3,1]  -  Position of the S/C in inertial frame [km]
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
%  mu      [1,1]  -  Gravitational constant of the observer [km^3/s^2]
%
%  R       [1,1]  -  Radius of the observer body [km]
%
%  n       [1,1]  - Truncation order of the zonal harmonics
%
%  JN      [n,1]  - Vector containing the first n zonal harmonics
%
%  et      [1,1]  - Initial time in seconds past J2000 
%
% OUTPUT ARGUMENTS:
%  a3B        [3x1]:   Perturbing acceleration in body-centered interial 
%                      frame [km/s]
% ------------------------------------------------------------------------
% CONTRIBUTOS: 
%  Alessio Derobertis
% ------------------------------------------------------------------------
% CHANGELOG:
%  04/04/2022 - First draft - ALessio Derobertis
%  11/04/2022 - Revision and bug fixing - Alessio Derobertis
%  10/01/2024 - New version for first release (Alessio Derobertis)
% ------------------------------------------------------------------------

if n > 1

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

    r_body = ROT*r_iner;

    % unpack the state
    rX = r_body(1);
    rY = r_body(2);
    rZ = r_body(3);
    rNorm = norm(r_body);

end

switch(n)

    case 0
       
        aGF_iner = [0,0,0];

    case 1

        aGF_iner = [0,0,0];

    case 2

        % J2 effect
        J2 = JN(1);

        aXj2 = - ((3*J2*mu*R^2*rX)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aYj2 = - ((3*J2*mu*R^2*rY)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aZj2 = - ((3*J2*mu*R^2*rZ)/(2*rNorm^5))*(3 - 5*(rZ^2/rNorm^2));

        aj2 = [aXj2; aYj2; aZj2];

        % total acceleration
        aGF_body = aj2;

        % convert to inertial ref. frame
        aGF_iner = ROT'*aGF_body;

    case 3

        % J2 effect
        J2 = JN(1);

        aXj2 = - ((3*J2*mu*R^2*rX)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aYj2 = - ((3*J2*mu*R^2*rY)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aZj2 = - ((3*J2*mu*R^2*rZ)/(2*rNorm^5))*(3 - 5*(rZ^2/rNorm^2));

        aj2 = [aXj2; aYj2; aZj2];

        % J3 effect
        J3 = JN(2);

        aXj3 = - ((5*J3*mu*R^3*rX)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aYj3 = - ((5*J3*mu*R^3*rY)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aZj3 = - ((5*J3*mu*R^3)/(2*rNorm^7))*(6*rZ^2 - 7*(rZ^4/rNorm^2) - (3/5)*rNorm^2);

        aj3 = [aXj3; aYj3; aZj3];

        % total acceleration
        aGF_body = aj2 + aj3;

        % convert to inertial ref. frame
        aGF_iner = ROT'*aGF_body;

    case 4

        % J2 effect
        J2 = JN(1);

        aXj2 = - ((3*J2*mu*R^2*rX)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aYj2 = - ((3*J2*mu*R^2*rY)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aZj2 = - ((3*J2*mu*R^2*rZ)/(2*rNorm^5))*(3 - 5*(rZ^2/rNorm^2));

        aj2 = [aXj2; aYj2; aZj2];

        % J3 effect
        J3 = JN(2);

        aXj3 = - ((5*J3*mu*R^3*rX)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aYj3 = - ((5*J3*mu*R^3*rY)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aZj3 = - ((5*J3*mu*R^3)/(2*rNorm^7))*(6*rZ^2 - 7*(rZ^4/rNorm^2) - (3/5)*rNorm^2);

        aj3 = [aXj3; aYj3; aZj3];

        % J4 effect
        J4 = JN(3);

        aXj4 = ((15*J4*mu*R^4*rX)/(8*rNorm^7))*(1 - (14*rZ^2)/rNorm^2 + (21*rZ^4)/rNorm^4);
        aYj4 = ((15*J4*mu*R^4*rY)/(8*rNorm^7))*(1 - (14*rZ^2)/rNorm^2 + (21*rZ^4)/rNorm^4);
        aZj4 = ((15*J4*mu*R^4*rZ)/(8*rNorm^7))*(5 - (70*rZ^2)/(3*rNorm^2) + (21*rZ^4)/rNorm^4);

        aj4 = [aXj4; aYj4; aZj4];

        % total acceleration
        aGF_body = aj2 + aj3 + aj4;
        
        % convert to inertial ref. frame
        aGF_iner = ROT'*aGF_body;

    case 5

        % J2 effect
        J2 = JN(1);

        aXj2 = - ((3*J2*mu*R^2*rX)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aYj2 = - ((3*J2*mu*R^2*rY)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aZj2 = - ((3*J2*mu*R^2*rZ)/(2*rNorm^5))*(3 - 5*(rZ^2/rNorm^2));

        aj2 = [aXj2; aYj2; aZj2];

        % J3 effect
        J3 = JN(2);

        aXj3 = - ((5*J3*mu*R^3*rX)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aYj3 = - ((5*J3*mu*R^3*rY)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aZj3 = - ((5*J3*mu*R^3)/(2*rNorm^7))*(6*rZ^2 - 7*(rZ^4/rNorm^2) - (3/5)*rNorm^2);

        aj3 = [aXj3; aYj3; aZj3];

        % J4 effect
        J4 = JN(3);

        aXj4 = ((15*J4*mu*R^4*rX)/(8*rNorm^7))*(1 - (14*rZ^2)/rNorm^2 + (21*rZ^4)/rNorm^4);
        aYj4 = ((15*J4*mu*R^4*rY)/(8*rNorm^7))*(1 - (14*rZ^2)/rNorm^2 + (21*rZ^4)/rNorm^4);
        aZj4 = ((15*J4*mu*R^4*rZ)/(8*rNorm^7))*(5 - (70*rZ^2)/(3*rNorm^2) + (21*rZ^4)/rNorm^4);

        aj4 = [aXj4; aYj4; aZj4];

        % J5 effect
        J5 = JN(4);

        aXj5 = ((3*J5*mu*R^5*rX*rZ)/(8*rNorm^9))*(35 - (210*rZ^2)/rNorm^2 + (231*rZ^4)/rNorm^4);
        aYj5 = ((3*J5*mu*R^5*rY*rZ)/(8*rNorm^9))*(35 - (210*rZ^2)/rNorm^2 + (231*rZ^4)/rNorm^4);
        aZj5 = ((3*J5*mu*R^5*rZ^2)/(8*rNorm^9))*(105 - (315*rZ^2)/rNorm^2 + (231*rZ^4)/rNorm^4) - (15*J5*mu*R^5)/(8*rNorm^7);

        aj5 = [aXj5; aYj5; aZj5];

        % total acceleration
        aGF_body = aj2 + aj3 + aj4 + aj5;

        % convert to inertial ref. frame
        aGF_iner = ROT'*aGF_body;

    case 6

        % J2 effect
        J2 = JN(1);

        aXj2 = - ((3*J2*mu*R^2*rX)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aYj2 = - ((3*J2*mu*R^2*rY)/(2*rNorm^5))*(1 - 5*(rZ^2/rNorm^2));
        aZj2 = - ((3*J2*mu*R^2*rZ)/(2*rNorm^5))*(3 - 5*(rZ^2/rNorm^2));

        aj2 = [aXj2; aYj2; aZj2];

        % J3 effect
        J3 = JN(2);

        aXj3 = - ((5*J3*mu*R^3*rX)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aYj3 = - ((5*J3*mu*R^3*rY)/(2*rNorm^7))*(3*rZ - 7*(rZ^3/rNorm^2));
        aZj3 = - ((5*J3*mu*R^3)/(2*rNorm^7))*(6*rZ^2 - 7*(rZ^4/rNorm^2) - (3/5)*rNorm^2);

        aj3 = [aXj3; aYj3; aZj3];

        % J4 effect
        J4 = JN(3);

        aXj4 = ((15*J4*mu*R^4*rX)/(8*rNorm^7))*(1 - (14*rZ^2)/rNorm^2 + (21*rZ^4)/rNorm^4);
        aYj4 = ((15*J4*mu*R^4*rY)/(8*rNorm^7))*(1 - (14*rZ^2)/rNorm^2 + (21*rZ^4)/rNorm^4);
        aZj4 = ((15*J4*mu*R^4*rZ)/(8*rNorm^7))*(5 - (70*rZ^2)/(3*rNorm^2) + (21*rZ^4)/rNorm^4);

        aj4 = [aXj4; aYj4; aZj4];

        % J5 effect
        J5 = JN(4);

        aXj5 = ((3*J5*mu*R^5*rX*rZ)/(8*rNorm^9))*(35 - (210*rZ^2)/rNorm^2 + (231*rZ^4)/rNorm^4);
        aYj5 = ((3*J5*mu*R^5*rY*rZ)/(8*rNorm^9))*(35 - (210*rZ^2)/rNorm^2 + (231*rZ^4)/rNorm^4);
        aZj5 = ((3*J5*mu*R^5*rZ^2)/(8*rNorm^9))*(105 - (315*rZ^2)/rNorm^2 + (231*rZ^4)/rNorm^4) - (15*J5*mu*R^5)/(8*rNorm^7);

        aj5 = [aXj5; aYj5; aZj5];

        % J6 effect
        J6 = JN(5);

        aXj6 = - ((J6*mu*R^6*rX)/(16*rNorm^9))*(35 - (945*rZ^2)/rNorm^2 + (3465*rZ^4)/rNorm^4 - (3003*rZ^6)/rNorm^6);
        aYj6 = - ((J6*mu*R^6*rY)/(16*rNorm^9))*(35 - (945*rZ^2)/rNorm^2 + (3465*rZ^4)/rNorm^4 - (3003*rZ^6)/rNorm^6);
        aZj6 = - ((J6*mu*R^6*rZ)/(16*rNorm^9))*(245 - (2205*rZ^2)/rNorm^2 + (4851*rZ^4)/rNorm^4 - (3003*rZ^6)/rNorm^6);

        aj6 = [aXj6; aYj6; aZj6];

        % total acceleration
        aGF_body = aj2 + aj3 + aj4 + aj5 + aj6;

        % convert to inertial ref. frame
        aGF_iner = ROT'*aGF_body;

end


end