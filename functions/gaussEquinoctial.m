function dequi = gaussEquinoctial(equi, mu, acc_rtn)
% ---------------------------------------------------------------------
% DESCRIPTION:
% GAUSSEQUINOCTIAL - Compute derivative of keplerian elements as Gauss
% Equinoctial elements
% ---------------------------------------------------------------------
% Compute derivatives of keplerian elements in the form of Gauss
% Equinoctial elements. This avoids singularities for near-circular
% and near-equatorial orbits.
% ---------------------------------------------------------------------
% INPUTS:
% equi         : [6,1] - equinoctial elements of the satellite ordered as:
%                        [p, f, g, h, k, L]
% mu           : [1,1] - gravitational parameter of the attractor
%                        [km^3/s^2]
% acc_rtn      : [3,1] - distrubing acceleration in rtn frame [km/s^2]
% ---------------------------------------------------------------------
% OUTPUTS:
% dequi        : [6,1] - derivatives of the equinoctial elements
%                        ordered as: [dp, df, dg, dh, dk, dL]
% ---------------------------------------------------------------------
% CONTRIBUTORS: Alessio Derobertis
% ---------------------------------------------------------------------

% unpack equinoctial elements
p = equi(1);
f = equi(2);
g = equi(3);
h = equi(4);
k = equi(5);
L = equi(6);

% unpack the acceleration in rtn
R = acc_rtn(1);
T = acc_rtn(2);
N = acc_rtn(3);

% compute derivatives
sL = sin(L);
cL = cos(L);
sqrtpm = sqrt(p/mu);
s2 = 1 + h^2 + k^2;
w = 1 + f*cL + g*sL;
s2no2w = s2*N / (2*w);

dp = (2*p*T/w)*sqrtpm;
df = sqrtpm * (R*sL + ((w+1)*cL+f)*T/w - g*N*(h*sL-k*cL)/w);
dg = sqrtpm * (-R*cL + ((w+1)*sL+g)*T/w + g*N*(h*sL-k*cL)/w);
dh = sqrtpm * s2no2w * cL;
dk = sqrtpm * s2no2w * sL;
dL = sqrt(mu*p)*(w/p)^2 + sqrtpm*((h*sL-k*cL)*N)/w;

dequi = [dp; df; dg; dh; dk; dL];

end