function Pnm_norm = expansionFunc(n)
%--------------------------------------------------------------------------
% EXPANSIONFUNC - Create symbolic functions to be used in the harmonic
% expansion series
%--------------------------------------------------------------------------
% Compute Normalized Associated Legendre Polynomial from order 0 up to 
% order n+1 using symbolic manipulation, the normalization is a full 
% normalization based on the formula:
% Pnm_norm = P_nm/[(n-m)!*(2*n+1)*(2-delta_0m)/(n+m)!]^(1/2)
% where delta_0m is = 1 if m=0
%                   = 0 if m!=0
% Also compute cos(j*x), sin(j*x) and (Re/r)^(1+1) functions.
%--------------------------------------------------------------------------
% INPUTS:
% n            : [1,1] - Truncation degree for the legendre polynomials [-]
% -------------------------------------------------------------------------
% OUTPUTS:
% Pnm_norm     : [(n+2)*(n+3)/2,1] - vector of normalised associated 
%                                    legendre polynomial up to order n+1 
%                                    ordered as column:
%                                      Pnm_norm(0) = poly(n=0,m=0)
%                                      Pnm_norm(1) = poly(n=1,m=0)
%                                      Pnm_norm(2) = poly(n=1,m=1)
%                                      Pnm_norm(3) = poly(n=2,m=0)
%                                      ...
% -------------------------------------------------------------------------
% SOURCES: Konopliv, A. S., et al. (2013), The JPL lunar gravity field to 
%          spherical harmonic degree 660 from the GRAIL Primary Mission, 
%          J. Geophys. Res. Planets, 118, 1415–1434, doi:10.1002/jgre.20097
% -------------------------------------------------------------------------
% CONTRIBUTORS: Alessio Derobertis
% -------------------------------------------------------------------------
% CHANGELOG:   2023/09/26 - V1: First draft
%              2023/10/17 - V2: Added computation of other functions
%              2024/01/17 - V2: Matlab release
% -------------------------------------------------------------------------

% Rise the order by 1
n = n + 1;

% Define the number of polynomials
g = (n + 1) * (n + 2) / 2;
g = int32(g);

% Initialize the function vectors
x = sym('x');
Pnm_norm = sym(zeros(g, 1));

% Define polynomials
c = 1;
for i = 0:n
    for j = 0:i
        if j == 0
            k = 1;
        else
            k = 2;
        end
        leg = legendreP(i,x);
        legDiff = diff(leg,x,j);
        f = sqrt(factorial(i - j) * ((2 * i) + 1) * k / factorial(i + j));
        Pnm_norm(c) = ((-1)^j)*((1 - x^2)^(j/2))*legDiff*f;
        c = c + 1;
    end
end

% Convert symbolic expressions to MATLAB functions
Pnm_norm = matlabFunction(Pnm_norm, 'Vars', x);

end