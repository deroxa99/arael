function [t,y] = arael(init_cond,ref_sys,perturb,settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARAEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Accurate tRAjectory Evaluation tooL %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% AUTHOR: Alessio Derobertis
%-------------------------------------------------------------------------%
% GITHUB: 
%-------------------------------------------------------------------------%
% RELEASE NOTES: Arael is an orbit propagater that allows the user to     
% choose between different propatgation methods, making use of various 
% levels of approximation of the dynamics in order to allow for           
% high-fidelity or high-performance simulations.                          
%                                                                         
% The selectable methods are:
%  - 'hifi'
%  - 'approx'
%  - 'averaged'
%  - 'full'
%
%------------------------------- HIFI ------------------------------------%
%  - GRAVITY: implement a complete harmonic expansion for the gravity 
%    field in order to simulate the dynamics in the most advanced way 
%    possible. The integration is carried out using Gauss Equinoctial 
%    elements and the supported central bodies are:
%    - 'EARTH': expansion up to order 360
%    - 'MOON': expansion up to order 1200
%    To speed up the computation, pre-computed polynomials are available 
%    for expansions of order n = [25,50,60,75,100,150,200,300].
%
%  - THIRD BODY: Perturbation due to third bodies can be implemented for
%    the following bodies:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%  
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%------------------------------ APPROX -----------------------------------%
%  - GRAVITY: implement an approximated model for the gravity field, to 
%    speed up the simulation. The integration is carried out using Gauss 
%    Equinoctial elements and up to now it is possible to simulate only 
%    zonal harmonics effect for the
%    following central bodies:
%     - 'EARTH': zonal harmonics up to J6;
%     - 'MOON': zonal harmonics up to J2;
%
%  - THIRD BODY: Perturbation due to third bodies can be implemented for
%    the following bodies:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%----------------------------- AVERAGED ----------------------------------%
%  - GRAVITY: implement an averaged zonal harmonic model that allows to cut
%    out fast evolving dynamics, drastically reducing the integration time.
%    Integration is carried out in Keplerian Elements and supported central
%    bodies are:
%     - 'EARTH': zonal harmonics up to J4;
%
%  - THIRD BODY: Perturbation due to third bodies are also averaged to cut
%    out fast evolving dynamics. Up to now the only perturbing bodies that
%    can be taken into account are:
%     - 'MOON', 'SUN'
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%------------------------------- FULL ------------------------------------%
%  - GRAVITY: implement an N-body problem where all the attractor are
%    treated as point masses. The integration is carried out using the
%    cartesian state. In this case perturb.n is not used.
%
%  - THIRD BODY: In this mode, the list of third bodies given in input is
%    used to compute the N-body dynamics. Allowed bodies are:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%-------------------------------------------------------------------------%
% THIS FUNCTION MAKE USE OF THE SPICE TOOLKIT DEVELOPED BY NASA, LEARN
% MORE ABOUT SPICE AT: https://naif.jpl.nasa.gov/naif/toolkit.html
% To work properly, locate the 'kernels' folder in the working folder.
%-------------------------------------------------------------------------%
% INPUT:
% init_cond  [-] - Struc containintg:
%                  .x0:     [6,1] - Initial state in inertial reference
%                                   frame [km,km/s]
%                  .et:     [1,1] - Initial time in seconds past J2000
%                  .tSpan:  [N,1] - Integration time-span [s]
%  
% ref_sys    [-] - Struct containing:
%                  .inertial [-] - String that defines the reference 
%                                  system as one of the following:
%                                   - Sun: 'ECLIPJ2000'
%                                   - Earth: 'J2000'
%                                   - Moon: 'MOON_PA_INERTIAL', 
%                                           'MOON_ME_INERTIAL' equivalent
%                                           to MOON_PA and MOON_ME @ J2000
%                                           epoch (1-1-2000 12:00:00)
%                  .obs      [-] - String with the name of the central 
%                                  body. Supported ones are:
%                                   - 'SUN'
%                                   - 'EARTH'
%                                   - 'MOON'
%  
% perturb    [-] - Struct containing:
%                  .n  [1,1] - Truncation order of the gravity field. In
%                              case of zonal harmoncis only, will be 
%                              considered terms up to Jn. 
%                  .TB {n,-} - Cell array containing the names 
%                              of the third-bodies attractor. In case of
%                              the N-body problem, those are all the
%                              attractor except for the central body.
%  
% settings   [-] - Struct containing:
%                   .mode    [-]   - Select integration mode:
%                                     - 'hifi': high-fidelity gravity
%                                     - 'approx': approximated dynamics
%                                     - 'averaged': averaged dynamics
%                                     - 'full': n_body problem
%                   .rel_rol [1,1] - Relative tolerance [-]
%                   .abs_tol [1,1] - Absolute tolerane [-]
%
% OUTPUT:
%
% t         [N,1] - Time span of the solution [s]
% y         [N,6] - Integrated state [km,km/s]
% 
%
%-------------------------------------------------------------------------%
% CHANGELOG: 2024/10/11 - Official relese for Matlab (Alessio Derobertis)
%-------------------------------------------------------------------------%

%%% add subfolders
addpath(genpath("arael"))

%%% load kernels
cspice_furnsh('arael/utils/meta_arael.tm');

%%% integrator options
options = odeset('RelTol',settings.rel_tol,'AbsTol',settings.abs_tol);

%%% Allowed settings

% reference system
ref_sys_allowed = {'ECLIPJ2000';'J2000';'MOON_PA_INERTIAL';'MOON_ME_INERTIAL'};
obs_allowed = {'SUN';'EARTH';'MOON'};

% gravity
n_hifi = [360;    % Earth
    1200];        % Moon

n_approx = [6;    % Earth
       2];        % Moon

n_averaged = [4]; % Earth

% third bodies
tb_allowed = {'SUN';
    'MERCURY';
    'VENUS';
    'EARTH';
    'MOON';
    'MARS BARYCENTER';
    'JUPITER BARYCENTER';
    'SATURN BARYCENTER';
    'URANUS BARYCENTER';
    'NEPTUNE BARYCENTER';
    'PLUTO BARYCENTER'};


%%% check for anomalies

% ref. system
if ismember(ref_sys.obs,obs_allowed) ~= 1
    fprintf('\n')
    fprintf('ERROR: invalid observer body \n')
    return
end

if ismember(ref_sys.inertial,ref_sys_allowed) ~= 1
    fprintf('\n')
    fprintf('ERROR: invalid reference system \n')
    return
end

% third-bodies
for i = 1:length(perturb.TB)
    if ismember(perturb.TB{i},tb_allowed) ~= 1
        fprintf('\n')
        fprintf(strcat('ERROR: ', perturbation.TB(i), ' is an invalid third-body \n'))
        return
    end
end

if ismember(ref_sys.obs,perturb.TB) == 1
    fprintf('\n')
    fprintf('ERROR: the observer body cannot be included in the third-body perturbators \n')
    return
end


switch settings.mode

    case 'hifi'

    case 'approx'
        % check for allowed observers
        obs_allowed_approx = {'EARTH';'MOON'};

        if ismember(ref_sys.obs,obs_allowed_approx) ~= 1
            fprintf('\n')
            fprintf('ERROR: invalid observer body \n')
            return
        end

        % compute the gravitational constant and radius
        mu = cspice_bodvrd(ref_sys.obs,'GM',1);
        Rp = mean(cspice_bodvrd(ref_sys.obs,'RADII',3));

        % GRAVITY
        aZH = @(t,r) pertZH(t,r,ref_sys.inertial,ref_sys.obs,mu,Rp,perturb.n,init_cond.et);

        % THIRD-BODY
        aTB = @(t,r) pertTB(t,r,perturb.TB,init_cond.et,ref_sys.inertial,ref_sys.obs);

        % PERTURBING ACCELERATION
        aTOT = @(t,x) aZH(t,x(1:3)) + aTB(t,x(1:3));

        % INITIAL CONDITION
        equi0 = car2equi(init_cond.x0,mu);

        % integrate using Equinoctial elements
        [t,equi_s] = ode113(@(t,equi) approx_rhs(t,equi,aTOT,mu),init_cond.tSpan,equi0,options);
        
        % convert equinoctial elements into state
        y = zeros(length(t),6);
        for i = 1:length(t)
            y(i,:) = equi2car(equi_s(i,:),mu);
        end

    case 'averaged'

    case 'full'

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% APPROX
function dequi = approx_rhs(t,equi,acc_ijk,mu)

% convert state into cartesian elements
x = equi2car(equi,mu);

% convert acceleration in rtn
r_ijk = x(1:3);
v_ijk = x(4:6);
r_vers = r_ijk/norm(r_ijk);
rv_cross = cross(r_ijk, v_ijk);
n_vers = rv_cross/norm(rv_cross);
t_vers = cross(n_vers, r_vers);
IJK_2_RTN = [r_vers, t_vers, n_vers];
acc_rtn = IJK_2_RTN * acc_ijk(t,x);

% compute variation of equinoctial elements
dequi = gaussEquinoctial(equi,mu,acc_rtn);
end