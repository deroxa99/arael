function [aTB] = pertTB(t,rSC,TB,et,ref_sys,obs)
% ------------------------------------------------------------------------
% DESCRIPTION:
% PERT3B - Compute perturbing acceleration due to third body attraction
% ------------------------------------------------------------------------
% INPUT ARGUMENTS:
%  t       [1,1]  -  Integration time [s]
%
%  rSC     [3,1]  -  Position of the S/C in inertial frame [km] 
%
%  TB      [n,1]  -  Vector of string containing the names of the
%                    third-bodies attractor. Supported ones are:
%                     - 'SUN','MERCURY','VENUS','EARTH','MOON',
%                     'MARS BARYCENTER', JUPITER BARYCENTER',
%                     'SATURN BARYCENTER', 'URANUS BARYCENTER',
%                     'NEPTUNE BARYCENTER','PLUTO BARYCENTER'.
%
%  et      [1,1]  -  Initial time in seconds past J2000 
%
%  ref_sys [-]    -  String that defines the reference system as one of
%                    the following:
%                     - Earth: 'J2000'
%                     - Moon: 'MOON_PA', 'MOON_ME' @ J2000 
%                             (1-1-2000 24:00:00)
%
%  obs     [-]    -  String containing the name of the central body.
%                    Supported ones are:
%                     - 'EARTH'
%                     - 'MOON'
%
% OUTPUT ARGUMENTS:
%  a3B        [3x1]:   Perturbing acceleration in body-centered inertial 
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

%%% state of the 3B
nTB = length(TB);
rTB = zeros(3,nTB);
muTB = zeros(nTB,1);

for i = 1:nTB
    rTB(:,i) = cspice_spkpos(TB{i},t+et,ref_sys,'NONE',obs);
    muTB(i) = cspice_bodvrd(TB{i},'GM',1);
end

rSCNorm = vecnorm(rSC);

rTBNorm = vecnorm(rTB);

rRel = rTB - rSC;
rRelNorm = vecnorm(rRel);

aTB = [0;0;0];
Q = zeros(3,nTB);

for i = 1:nTB

    % aTB = aTB + muTB(i).*(rRel(:,i)/rRelNorm(i)^3 - rTB(:,i)/rTBNorm^3);

    %%% impelent Roy formulation for reduced error

    Q(:,i) = ((rSCNorm.^2 + 2*dot(rSC,rRel(:,i)))/(rTBNorm(i)^3*rRelNorm(i)^3))*((rTBNorm(i)^2 + rTBNorm(i)*rRelNorm(i) + rRelNorm(i)^2)./(rTBNorm(:,i) + rRelNorm(:,i)));

    aTB = aTB + muTB(i).*(rRel(:,i).*Q(:,i) - rSC./rTBNorm(i)^3);

end

end
