function [dkepTB] = avgTB(t,kep,TB,et,ref_sys,obs)
% ------------------------------------------------------------------------
% DESCRIPTION:
% AVG3B - Compute averaged perturbing acceleration due to third body 
% attraction
% ------------------------------------------------------------------------
% INPUT ARGUMENTS:
%  t       [1,1]  -  Integration time [s]
%
%  rSC     [3,1]  -  Position of the S/C in inertial frame [km] 
%
%  TB      [n,1]  -  Vector of string containing the names of the
%                    third-bodies attractor. Supported ones are:
%                     - 'EARTH'
%                     - 'MOON'
%                     - 'SUN'
%
%  et      [1,1]  -  Initial time in seconds past J2000 
%
%  ref_sys [-]    -  String that defines the reference system as one of
%                    the following:
%                     - Earth: 'J2000'
%                     - Moon: 'MOON_PA_INERTIAL', 'MOON_ME_INERTIAL'
%
%  obs     [-]    -  String containing the name of the central body.
%                    Supported ones are:
%                     - 'EARTH'
%                     - 'MOON'
%
% OUTPUT ARGUMENTS:
%  a3B        [3x1]:   Perturbing acceleration [km/s]
% ------------------------------------------------------------------------
% CONTRIBUTOS: 
%  Alessio Derobertis
% ------------------------------------------------------------------------
% CHANGELOG: 
%  04/04/2022 - First draft - ALessio Derobertis
%  11/04/2022 - Revision and bug fixing - Alessio Derobertis
%  10/01/2024 - New version for first release (Alessio Derobertis)
% ------------------------------------------------------------------------

%%% keplerian elements of the S/C
a = kep(1);
e = kep(2);
i = kep(3);
RAAN = kep(4);
omega = kep(5);
f0 = kep(6);

%%% state of the 3B
mu = cspice_bodvrd(obs,'GM',1);
nTB = length(TB);
sTB = zeros(6,nTB);
muTB = zeros(nTB,1);
massRatio = zeros(nTB,1);

for i = 1:nTB
    if strcmp(TB(i),'SUN')
        TB_temp = strcat(TB(i),' BARYCENTER');
        massRatio(i) = 1;
    else
        TB_temp = TB(i);
    end

    sTB(:,i) = cpsice_spkpos(TB_temp,t + et,ref_sys,'NONE',obs);
    muTB(i) = cspice_bodvrd(TB(i),'GM',1);
    massRatio(i) = muTB(i)/mu;
end


for j = 1:num3B

    kepTB = car2kep(sTB(1:3,j),v3B(4:6,j),mu);
    rTBnorm = norm(sTBB(1:3),i);

    aTB = kepTB(1);
    eTB = kepTB(2);
    iTB = kepTB(3);
    RAANTB = kepTB(4);
    omegaTB = kepTB(5);
    f0TB = kepTB(6);

    u = omega + f0;
    u3B = omegaTB + f0TB;

    deltaRAAN = RAAN - RAANTB;

    nP = sqrt(muP/a^3);
    n3B = sqrt(muAtt(j)/aTB^3);

    gamma = (n3B^2/nP)*(aTB/rTBnorm)^3*massRatio(j);
    s = sqrt(1 - e^2);

    %%% series expansion

    ApB = (1/4)*(1 + cos(i)^2)*(1 + cos(iTB)^2) + (1/2)*sin(i)^2*sin(iTB)^2 + ...
        (1/4)*(sin(2*i)*sin(2*iTB)*cos(deltaRAAN) + sin(i)^2*sin(iTB)^2*cos(2*deltaRAAN)) ...
        + (1/2)*(1 - (3/2)*sin(i)^2)*sin(iTB)^2*cos(2*u) ...
        + (1/8)*sin(i)^2*((1 + cos(i))^2*cos(2*(u - deltaRAAN)) + (1 - cos(i))^2*cos(2*(u + deltaRAAN))) + ...
        + (1/4)*sin(2*i)*sin(iTB)*((1 - cos(iTB))*cos(2*u + deltaRAAN) - (1 + cos(iTB))*cos(2*u - deltaRAAN));

    AmB = (1/2)*sin(i)^2*(1 - (3/2)*sin(iTB)^2)*cos(2*omega) + (1/4)*sin(2*iTB)*sin(i)*((1 - cos(i))*cos(2*omega - deltaRAAN) ...
        - (1 + cos(i))*cos(2*omega + deltaRAAN)) + (1/8)*sin(iTB)^2*((1 - cos(i))^2*cos(2*(omega - deltaRAAN)) + (1 + cos(i))^2*cos(2*(omega + deltaRAAN))) ...
        + (3/8)*sin(i)^2*sin(iTB)^2*(cos(2*(u - omega)) + cos(2*(u + omega))) ...
        + (1/4)*sin(iTB/2)^4*(sin(i/2)^4*cos(2*(u - omega + deltaRAAN)) + cos(i/2)^4*cos(2*(u + omega + deltaRAAN))) ...
        + (1/4)*cos(iTB/2)^4*(sin(i/2)^4*cos(2*(u + omega - deltaRAAN)) + cos(i/2)^4*cos(2*(u - omega - deltaRAAN))) ...
        + sin(i)*sin(iTB)*sin(iTB/2)^2*(sin(i/2)^2*cos(2*u - 2*omega + deltaRAAN) - cos(i/2)^2*cos(2*u + 2*omega + deltaRAAN)) ...
        + sin(i)*sin(iTB)*cos(iTB/2)^2*(cos(i/2)^2*cos(2*u - 2*omega - deltaRAAN) - sin(i/2)^2*cos(2*u + 2*omega - deltaRAAN));

    IPO = -(1/2)*sin(i)*((cos(i)*sin(2*iTB)*sin(deltaRAAN) + sin(i)*sin(iTB)^2*sin(2*deltaRAAN)) ...
        - (1/2)*sin(i)*((1 + cos(iTB))^2*sin(2*(u3B - deltaRAAN)) - (1 - cos(iTB))^2*sin(2*(u3B + deltaRAAN))) ...
        + cos(i)*sin(iTB)*((1 + cos(iTB))*sin(2*u3B - deltaRAAN) + (1 - cos(iTB)*sin(2*u3B + deltaRAAN))));

    IMS = -sin(i)*sin(iTB)^2*cos(2*omega) + (1/4)*sin(2*iTB)*((1 - cos(i))*cos(2*omega - deltaRAAN) ...
        - (1 + cos(i))*cos(2*omega + deltaRAAN)) ...
        + (1/2)*sin(i)*sin(iTB)^2*(cos(2*(u3B - omega)) + cos(2*(u3B + omega))) ...
        + (1/4)*sin(iTB)*((1 + cos(i))*(1 + cos(iTB))*cos(2*u3B - 2*omega - deltaRAAN)) ...
        + (1 - cos(i))*(1 - cos(iTB))*cos(2*u3B - 2*omega + deltaRAAN)  ...
        - (1 - cos(i))*(1 + cos(iTB))*cos(2*u3B + 2*omega - deltaRAAN)  ...
        - (1 + cos(i))*(1 - cos(iTB))*cos(2*u3B + 2*omega + deltaRAAN);

    IMO = (1/4)*(sin(iTB)^2*((1 - cos(i))^2*sin(2*(omega - deltaRAAN) - (1 + cos(i))^2*sin(2*(omega + deltaRAAN)))) ...
        + sin(i)*sin(2*iTB)*((1 + cos(i))*sin(2*omega + deltaRAAN) + (1 - cos(i))*sin(2*omega - deltaRAAN)) ...
        + (1/2)*((1 + cos(i))^2*(1 + cos(iTB))^2*sin(2*(u3B - omega - deltaRAAN)) ...
        - (1 - cos(i))^2*(1 - cos(iTB))^2*sin(2*(u3B - omega + deltaRAAN)) ...
        + (1 - cos(i))^2*(1 + cos(iTB))^2*sin(2*(u3B + omega - deltaRAAN)) ...
        - (1 + cos(i))^2*(1 - cos(iTB))^2*sin(2*(u3B + omega + deltaRAAN))) ...
        + sin(i)*sin(iTB)*((1 + cos(i))*(1 + cos(iTB))*sin(2*u3B - 2*omega - deltaRAAN) ...
        - (1 - cos(i))*(1 - cos(iTB))*sin(2*u3B - 2*omega + deltaRAAN) ...
        - (1 - cos(i))*(1 + cos(iTB))*sin(2*u3B + 2*omega - deltaRAAN) ...
        + (1 + cos(i))*(1 - cos(iTB))*sin(2*u3B + 2*omega + deltaRAAN)));

    IPS = sin(i)*sin(iTB)^2 + (1/2)*cos(i)*sin(2*iTB)*cos(deltaRAAN) - sin(i)*sin(iTB)^2*cos(2*u3B) ...
        + (1/2)*cos(i)*sin(iTB)*((1 - cos(iTB))*cos(2*u3B + deltaRAAN) - (1 + cos(iTB))*cos(2*u3B - deltaRAAN));

    IPC = cos(i)*(1 - (1/2)*sin(iTB)^2) + (1/2)*sin(i)*sin(2*iTB)*cos(deltaRAAN) ...
        - (1/2)*cos(i)*sin(iTB)^2*cos(2*deltaRAAN) + (1/2)*cos(i)*sin(iTB)^2*cos(2*u3B) ...
        + (1/2)*sin(i)*sin(iTB)*((1 - cos(iTB))*cos(2*u3B + deltaRAAN) ...
        - (1 + cos(iTB))*cos(2*u3B - deltaRAAN)) - (1/4)*cos(i)*((1 + cos(iTB))^2*cos(2*(u3B - deltaRAAN)) + ...
        (1 - cos(iTB))^2*cos(2*(u3B + deltaRAAN)));

    IMW = (1/2)*sin(i)^2*(3*sin(iTB)^2 - 2)*sin(2*omega) ...
        - (1/4)*sin(iTB)^2*((1 - cos(i))^2*sin(2*(omega - deltaRAAN)) + (1 + cos(i))^2*sin(2*(omega + deltaRAAN))) ...
        + (1/2)*sin(i)*sin(2*iTB)*((1 + cos(i))*sin(2*omega + deltaRAAN) - (1 - cos(i))*sin(2*omega - deltaRAAN)) ...
        + (1/8)*((1 + cos(i))^2*(1 + cos(iTB))^2*sin(2*(u3B - omega - deltaRAAN)) ...
        + (1 - cos(i))^2*(1 - cos(iTB))^2*sin(2*(u3B - omega + deltaRAAN))  ...
        - (1 - cos(i))^2*(1 + cos(iTB))^2*sin(2*(u3B + omega - deltaRAAN))  ...
        - (1 + cos(i))^2*(1 - cos(iTB))^2*sin(2*(u3B + omega + deltaRAAN)))  ...
        + (3/4)*sin(i)^2*sin(iTB)^2*(sin(2*(u3B - omega)) - sin(2*(u3B + omega))) ...
        + (1/2)*sin(i)*sin(iTB)*((1 + cos(i))*(1 + cos(iTB))*sin(2*u3B - 2*omega - deltaRAAN) ...
        + (1 - cos(i))*(1 - cos(iTB))*sin(2*u3B - 2*omega + deltaRAAN)  ...
        + (1 - cos(i))*(1 + cos(iTB))*sin(2*u3B + 2*omega - deltaRAAN)  ...
        + (1 + cos(i))*(1 - cos(iTB))*sin(2*u3B + 2*omega + deltaRAAN));

    IMC = -(1/4)*sin(i)*sin(2*iTB)*(cos(2*omega - deltaRAAN) + cos(2*omega + deltaRAAN)) ...
        + cos(i)*((1/2)*sin(iTB)^2 - 1)*cos(2*omega) ...
        + (1/4)*sin(iTB)^2*((1 + cos(i))*cos(2*(omega + deltaRAAN)) - (1 - cos(i))*cos(2*(omega - deltaRAAN))) ...
        + (1/8)*((1 + cos(i))*(1 + cos(iTB))^2*cos(2*(u3B - omega - deltaRAAN)) ...
        - (1 - cos(i))*(1 - cos(iTB))^2*cos(2*(u3B - omega + deltaRAAN))...
        - (1 - cos(i))*(1 + cos(iTB))^2*cos(2*(u3B + omega - deltaRAAN))...
        + (1 + cos(i))*(1 - cos(iTB))^2*cos(2*(u3B + omega + deltaRAAN))) ...
        + (1/4)*sin(i)*sin(iTB)*(1 + cos(iTB))*(cos(2*u3B - 2*omega - deltaRAAN) + cos(2*u3B + 2*omega - deltaRAAN)) ...
        - (1/4)*sin(i)*sin(iTB)*(1 - cos(iTB))*(cos(2*u3B - 2*omega + deltaRAAN) + cos(2*u3B + 2*omega + deltaRAAN)) ...
        - (1/4)*cos(i)*sin(iTB)^2*(cos(2*(u3B - omega)) + cos(2*(u3B + omega)));

    %%% variation of keplerian elements

    da = 0;

    de = - (15/8)*e*gamma*s*IMW;

    di = - (3/4)*((gamma*csc(i))/s)*((1 + 3/2*e^2)*IPO + (5/2)*e^2*(IMO - cos(i)*IMW));

    dRAAN = (3/4)*(gamma/s)*((1 + (3/2)*e^2)*((cos(i)/sin(i))*IPS - IPC) + (5/2)*e^2*((cos(i)/sin(i))*IMS - IMC));

    domega = (3/2)*s*gamma*((3/2)*(ApB) - 1 + (5/2)*AmB) - cos(i)*dRAAN;

    df0 = - 2*gamma*((1 + (3/2)*e^2)*((3/2)*ApB - 1) + (15/4)*e^2*AmB) - (3/2)*s^2*gamma*(3/2*ApB - 1 + (5/2)*AmB);

    dkepTB = dkepTB + [da; de; di; dRAAN; domega; df0];

end

end
