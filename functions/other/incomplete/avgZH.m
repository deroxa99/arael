function [dyZH] = avgZH(t,kepSC,J2,J3,J4,muP,RP,options)
% ----------------------------------------------------------------------
% DESCRIPTION:
% Compute averaged variations of Keplerian Elements form perturbation due to
% gravity zonal harmonics [J2,J3,J4] 
% ----------------------------------------------------------------------
% INPUT ARGUMENTS:
%   t             [1X1]   time                    
%   kepSC         [6X1]   Vector containing the keplerian elements of the
%                         spacecraft, in this order: 
%                         a - Semi-major Axis [km]
%                         e - Eccentricity [-]
%                         i - Inclination [rad]
%                         RAAN - Right Ascesion of Ascending Node [rad]
%                         omega - Pericentre Anomaly [rad]
%                         f0 - True Anomaly [rad]
%   J2            [1X1]   J2 parameter for zonal harmonics [-]
%   J3            [1X1]   J3 parameter for zonal harmonics [-]
%   muP           [1X1]   Planetary Gravitational constant of the 
%                         planet [km^3/s^2]
%   RP            [1X1]   Mean radius of the planet [km]
%   options       [1X1]   Number of harmonics to include (min 1, max 4)           
%
% OUTPUT ARGUMENTS:
%   dyZH          [6x1]  Derivative of keplerian elements
% ----------------------------------------------------------------------
% REFERENCES:
%  - Applied Orbit Perturbation and Maintenance - Chia-Chun "George" Chao - 2005
%  - A Semi-Analytical Solution for the Motion of a Low Altitude Earth
%    Satellite Under j2-gravity and air drag perturbations - H.A Embary,
%    A.H. Ibrahim, I.A Hassa M.N Ismail - 2021
% ----------------------------------------------------------------------
% CONTRIBUTOS: 
%   Alessio Derobertis
% ----------------------------------------------------------------------
% CHANGELOG:
%   04/04/2022 - First draft - ALessio Derobertis
%   11/04/2022 - Added J4 influence - Alessio Derobertis
%   01/10/2024 - First release for Matlab (Alessio Derobertis)
% ----------------------------------------------------------------------

%%% keplerian elements of the S/C

a = kepSC(1);
e = kepSC(2);
i = kepSC(3);
RAAN = kepSC(4);
omega = kepSC(5);
f0 = kepSC(6);

p = a*(1 - e^2);

n = sqrt(muP/a^3);

switch options
    case 1
        da = 0;

        de = 0;

        di = 0;

        dRAAN = 0;

        domega = 0;

        df0 = n;

        dyZH = [da; de; di; dRAAN; domega; df0];

    case 2
        da = 0;

        de = 0;

        di = 0;

        dRAAN = -(3/2)*n*J2*(RP/p)^2*cos(i);

        domega = (3/4)*n*J2*(RP/p)^2*(5*cos(i)^2 - 1);

        df0 = n + (3/(4*(1 - e^2)^(3/2)))*n*J2*(RP/p)^2*(3*cos(i)^2 - 1);

        dyZH = [da; de; di; dRAAN; domega; df0];

    case 3
        %%% variation of keplerian elements (Chia Chung form)

        % da = 0;
        %
        % de = -(3/8)*n*J3*(RP/p)^3*sin(i)*(4 - 5*sin(i)^2)*(1 - e^2)*cos(omega);
        %
        % di = (3/8)*n*J3*(RP/p)^3*cos(i)*(4 - 5*sin(i)^2)*e*cos(omega);
        %
        % dRAAN = -(3/2)*n*J2*(RP/p)^2*cos(i) ...
        %    - (3/8)*n*J3*(RP/p)^3*(15*sin(i)^2 - 4)*(e*cot(i)*sin(omega));
        %
        % domega = (3/4)*n*J2*(RP/p)^2*(4 - 5*sin(i)^2) ...
        %    + (3/8)*n*J3*(RP/p)^3*(((4 - 5*sin(i)^2)*(sin(i)^2 - e^2*cos(i)^2))/(e*sin(i)) + 2*sin(i)*(13 - 15*sin(i)^2)*e)*sin(omega);
        %
        % df0 = n*(1 + (3/2)*J2*(RP/p)^2*(1 - (3/2)*sin(i)^2)*(1 - e^2)^0.5) - ...
        %    (3/8)*n*J3*(RP/p)^3*sin(i)*(4 - 5*sin(i)^2)*(1 - 4*e^2)*(((1 - e^2)^0.5)/e)*sin(omega);
        %
        % dyZH = [da; de; di; dRAAN; domega; df0];

        %%% variation of keplerian elements (Embary form)

        da = 0;

        de = -(3/16)*n*J3*(RP/p)^3*(3 + 5*cos(2*i))*(1 - e^2)*cos(omega)*sin(i);

        di = - (3/16)*n*e*J3*(RP/p)^3*(3 + 5*cos(2*i))*cos(i)*cos(omega);

        dRAAN = -(3/2)*n*J2*(RP/p)^2*cos(i) ...
            - (3/15)*n*J3*(RP/p)^3*(cos(i) + 15*cos(3*i))*csc(i)*sin(omega);

        domega = (3/4)*n*J2*(RP/p)^2*(5*cos(i)^2 - 1) ...
            + (3/(e*64))*n*J3*(RP/p)^3*sin(omega)*csc(i)*(-1 - 3*e^2 - 4*cos(2*i) + 5*(1 + 7*e^2)*cos(4*i));

        df0 = n + (3/(4*(1 - e^2)^(3/2)))*n*J2*(RP/p)^2*(3*cos(i)^2 - 1) ...
            + (3/32)*n*J3*(RP/a)^3*sin(omega)*(((1 + 12*e^2)*(sin(i) + 5*sin(3*i)))/(e*(1 - e^2)^(5/2)));

        dyZH = [da; de; di; dRAAN; domega; df0];

    case 4
        da = 0;

        de = -(3/16)*n*J3*(RP/p)^3*(3 + 5*cos(2*i))*(1 - e^2)*cos(omega)*sin(i) ...
            + (15/64)*e*n*J4*(RP/p)^4*sin(i)^2*sin(2*omega)*(5 + 7*cos(2*i))*(1 - e^2);

        di = - (3/16)*n*e*J3*(RP/p)^3*(3 + 5*cos(2*i))*cos(i)*cos(omega) ...
            + (15/256)*e^2*n*J4*(RP/p)^4*sin(2*omega)*(10*sin(2*i) + 7*sin(4*i));

        dRAAN = -(3/2)*n*J2*(RP/p)^2*cos(i) ...
            - (3/15)*n*J3*(RP/p)^3*(cos(i) + 15*cos(3*i))*csc(i)*sin(omega) ...
            + (15/128)*n*J4*(RP/a)^4*((2 + 3*e^2)*(9*cos(i) + 7*cos(3*i)) -2*e^2*cos(2*omega)*(5*cos(i) + 7*cos(3*i)));

        domega = (3/4)*n*J2*(RP/p)^2*(5*cos(i)^2 - 1) ...
            + (3/(e*64))*n*J3*(RP/p)^3*sin(omega)*csc(i)*(-1 - 3*e^2 - 4*cos(2*i) + 5*(1 + 7*e^2)*cos(4*i)) ...
            + (15/1024)*n*J4*(RP/a)^4*(-27*(4 + 5*e^2) + 2*(-6 + 5*e^2)*cos(2*omega) + 4*cos(2*i)*(-52 - 63*e^2 + 2*(-2 + 7*e^2)*cos(2*omega)));

        df0 = n + (3/(4*(1 - e^2)^(3/2)))*n*J2*(RP/p)^2*(3*cos(i)^2 - 1) ...
            + (3/32)*n*J3*(RP/a)^3*sin(omega)*(((1 + 12*e^2)*(sin(i) + 5*sin(3*i)))/(e*(1 - e^2)^(5/2))) ...
            + (15/1024)*n*J4*(RP/a)^4*(((8 + 9*e^2)*(9 + 20*cos(2*i) + 35*cos(4*i)) + 8*(2 + 15*e^2)*(5 + 7*cos(2*i))*cos(2*omega)*sin(i)^2)/(1 - e^2)^(7/2));

        dyZH = [da; de; di; dRAAN; domega; df0];

end

end