function [t,y] = arael(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ARAEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Adaptable tRAjectory Evaluation tooL %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
% AUTHOR: Alessio Derobertis
%-------------------------------------------------------------------------%
% GITHUB: https://github.com/deroxa99/arael
%-------------------------------------------------------------------------%
% RELEASE NOTES: Arael is an orbit propagater that allows the user to     
% choose between different propatgation methods, making use of various 
% levels of approximation of the dynamics in order to allow for a trade off
% between high-fidelity or high-performance simulations.                        
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
%    - 'EARTH': expansion up to order 360 (EGM96)
%               - Source: https://cddis.nasa.gov/926/egm96/
%    - 'MOON': expansion up to order 1200 (GRGM1200A)
%               - Source: https://pgda.gsfc.nasa.gov/products/50
%    - 'MARS': expansion up to order 80 (GGM1025A)
%               - Source: https://pds-ppi.igpp.ucla.edu/search/view/
%                         ?f=yes&id=pds://PPI/MGS-M-RSS-5-SDP-V1.0/DATA/
%                         RS_SHA/GGM1025A&o=1
%    - 'VENUS': expansion up to order 180 (MGN180u)
%               - Source: https://pds-geosciences.wustl.edu/mgn/
%                         mgn-v-rss-5-gravity-l2-v1/mg_5201/gravity/
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
%  - SOLAR RADIATION PRESSURE: Perturbation due to solar radiation hitting
%    the surface of the satellite. Umbra and penumbra regions due to the
%    shadow of the primary attractor are also taken into account.
%
%------------------------------ APPROX -----------------------------------%
%  - GRAVITY: implement an approximated model for the gravity field, to 
%    speed up the simulation. The integration is carried out using Gauss 
%    Equinoctial elements and up to now it is possible to simulate only 
%    zonal harmonics effect for the
%    following central bodies:
%     - 'EARTH': zonal harmonics up to J6;
%     - 'MOON': zonal harmonics up to J4;
%     - 'MARS': zonal harmonics up to J4;
%     - 'VENUS': zonal harmonics up to J4;
%
%  - THIRD BODY: Perturbation due to third bodies can be implemented for
%    the following bodies:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: Perturbation due to solar radiation hitting
%    the surface of the satellite. Umbra and penumbra regions due to the
%    shadow of the primary attractor are also taken into account.
%
%----------------------------- AVERAGED ----------------------------------%
%  - GRAVITY: not yet implemented
%
%  - THIRD BODY: not yet implemented
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%------------------------------- FULL ------------------------------------%
%  - GRAVITY: implement an N-body problem where all the attractor are
%    treated as point masses. The integration is carried out using the
%    cartesian state. In this case perturb.n is not used. Allowed observer
%    bodies are:
%    - 'SUN'
%
%  - THIRD BODY: In this mode, the list of third bodies given in input is
%    used to compute the N-body dynamics. Allowed bodies are:
%    - 'SUN','MERCURY','VENUS','EARTH','MOON','MARS BARYCENTER',
%      'JUPITER BARYCENTER','SATURN BARYCENTER', 'URANUS BARYCENTER',
%      'NEPTUNE BARYCENTER','PLUTO' BARYCENTER.
%
%  - AIR DRAG: not yet implemented
%
%  - SOLAR RADIATION PRESSURE: Perturbation due to solar radiation hitting
%    the surface of the satellite. Shadow is not considered since the S/C
%    is supposed to orbit the Sun.
%
%-------------------------------------------------------------------------%
% THIS FUNCTION MAKE USE OF THE SPICE TOOLKIT DEVELOPED BY NASA, LEARN
% MORE ABOUT SPICE AT: https://naif.jpl.nasa.gov/naif/toolkit.html
%-------------------------------------------------------------------------%
% NOTE: input is CAPS SENSITIVE, no need for ordered input
%-------------------------------------------------------------------------%
% INPUT: 
% init_cond  [-] - Struc containintg:
%                  .x0:     [6,1] - Initial state in inertial reference
%                                   frame [km,km/s]
%                  .utc0:   [-]   - Initial time UTC, formatted as:
%                                   'yyyy-mm-dd hh:mm:ss.sss UTC'
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
%                                  observer body. Supported ones depends on
%                                  the mode (see each mode description).
%  
% perturb    [-] - Struct containing:
%                  .n  [1,1] - Truncation order of the gravity field. In
%                              case of zonal harmoncis only, will be 
%                              considered terms up to Jn. 
%                              Set as 0 or 1 for point mass.
%                              (default = 0).
%                  .TB {n,-} - Cell array containing the names 
%                              of the third-bodies attractor. In case of
%                              the N-body problem, those are all the
%                              attractor except for the central body.
%                              Leave empty if no third-body perturbation is
%                              taken into account (default = {}).
%                  .SRP [-]  - Set as:
%                               - 'on': SRP is included
%                               - 'off': SRP is not included
%                               (default = 'off').
%
% spacecraft [-] - Struct containing: 
%                  .m [1,1]  - Mass of the S/C [km] (default = 850 kg)
%                  .A [1,1]  - Cross-section of the S/C [m^2] 
%                              (default = 2 m^2)
%                  .cR [1,1] - Reflectivity coefficient of the S/C [-]
%                              (default = 1.8)
%  
% settings   [-] - Struct containing:
%                  .mode    [-]   - Select integration mode:
%                                   - 'hifi': high-fidelity gravity
%                                   - 'approx': approximated dynamics
%                                   - 'averaged': averaged dynamics
%                                   - 'full': n_body problem
%                                  (default is 'hifi')
%                  .rel_tol [1,1] - Relative tolerance (default = 1e-09)
%                  .abs_tol [1,1] - Absolute tolerane (default = 1e-09)
%
% OUTPUT:
%
% t          [N,1] - Time span of the solution [s]
% y          [N,6] - Integrated state [km,km/s] 
% 
%-------------------------------------------------------------------------%
% CHANGELOG: 2024/10/11 - Official relese for Matlab (Alessio Derobertis)
%-------------------------------------------------------------------------%
% TO DO: - debug MARS and VENUS using GMAT
%        - check zonal coefficients for MOON,MARS,VENUS
%        - add 'average'
%        - add DRAG
%-------------------------------------------------------------------------%

%%% add subfolders
addpath(genpath("arael"))

%%% load kernels
cspice_furnsh('arael/utils/meta_arael.tm');

%%% define input
possible_in = {'init_cond';'ref_sys';'perturb';'spacecraft';'settings'};

for i = 1:nargin
    if strcmp(inputname(i),possible_in{1}) == 1
        init_cond = varargin{i};
        break
    else 
        init_cond = struct;
    end
end

for i = 1:nargin
    if strcmp(inputname(i),possible_in{2}) == 1
        ref_sys = varargin{i};
        break
    else
        ref_sys = struct;
    end
end

for i = 1:nargin
    if strcmp(inputname(i),possible_in{3}) == 1
        perturb = varargin{i};
        break
    else
        perturb = struct;
    end
end

for i = 1:nargin
    if strcmp(inputname(i),possible_in{4}) == 1
        spacecraft = varargin{i};
        break
    else
        spacecraft = struct;
    end
end

for i = 1:nargin
    if strcmp(inputname(i),possible_in{5}) == 1
        settings = varargin{i};
        break
    else 
        settings = struct;
    end
end

%%% define error output
t = 0;
y = 0;

%%% Allowed settings

% reference system
ref_sys_allowed = {'ECLIPJ2000';'J2000';'MOON_PA_INERTIAL';'MOON_ME_INERTIAL';
    };
obs_allowed = {'SUN';'EARTH';'MOON';'MARS';'VENUS'};

% gravity
n_hifi = [360;    % Earth
    1200;         % Moon
    80;           % Mars
    180];         % Venus      

n_approx = [6;    % Earth
       4          % Moon
       4          % Mars
       4];        % Venus

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

% SRP
SRP_allowed = {'on';'off'};

%%% check for anomalies

% check input
if isfield(init_cond,'x0') ~= 1
    fprintf('ERROR: initial state is missing (init_cond.x0)\n')
    return
end

if isfield(init_cond,'utc0') ~= 1
    fprintf('ERROR: initial time is missing (init_cond.utc0)\n')
    return
end

if isfield(init_cond,'tSpan') ~= 1
    fprintf('ERROR: integration time span is missing (init_cond.tSpan)\n')
    return
end

if isfield(ref_sys,'inertial') ~= 1
    fprintf('ERROR: inertial reference system is missing (ref_sys.inertial)\n')
    return
end

if isfield(ref_sys,'obs') ~= 1
    fprintf('ERROR: observer body is missing (ref_sys.obs)\n')
    return
end

if isfield(perturb,'n') ~= 1
    fprintf('WARNING: perturb.n is missing, set to default\n')
    perturb.n = 0
end

if isfield(perturb,'TB') ~= 1
    fprintf('WARNING: perturb.TB is missing, set to default\n')
    perturb.TB = {}
end

if isfield(perturb,'SRP') ~= 1
    fprintf('WARNING: perturb.SRP is missing, set to default\n')
    perturb.SRP = 'off'
end

if isfield(spacecraft,'m') ~= 1
    fprintf('WARNING: spacecraft.m is missing, set to default\n')
    spacecraft.m = 850
end

if isfield(spacecraft,'A') ~= 1
    fprintf('WARNING: spacecraft.A is missing, set to default\n')
    spacecraft.A = 2
end

if isfield(spacecraft,'cR') ~= 1
    fprintf('WARNING: spacecraft.cR is missing, set to default\n')
    spacecraft.cR = 1.8
end

if isfield(settings,'mode') ~= 1
    fprintf('WARNING: settings.mode missing, set to default\n')  
    settings.mode = 'hifi'
end

if isfield(settings,'rel_tol') ~= 1
    fprintf('WARNING: settings.rel_tol is missing, set to default\n')
    settings.rel_tol = 1e-09
end

if isfield(settings,'abs_tol') ~= 1
    fprintf('WARNING: settings.abs_tol is missing, set to default\n')
    settings.abs_tol = 1e-09
end

%%% integrator options
options = odeset('RelTol',settings.rel_tol,'AbsTol',settings.abs_tol);

%%% initial time in et
et = cspice_str2et(init_cond.utc0);

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

if ismember(perturb.SRP,SRP_allowed) ~= 1
    fprintf('\n')
    fprintf('ERROR: invalid solar radiation pressure settings \n')
    return
end


switch settings.mode

    case 'hifi'
        % check for allowed observers
        obs_allowed_approx = {'EARTH';'MOON';'MARS';'VENUS'};

        if ismember(ref_sys.obs,obs_allowed_approx) ~= 1
            fprintf('\n')
            fprintf('ERROR: invalid observer body \n')
            return
        end

        % compute the gravitational constant and radius
        mu = cspice_bodvrd(ref_sys.obs,'GM',1);
        Rp = mean(cspice_bodvrd(ref_sys.obs,'RADII',3));
        
        muTB = zeros(length(perturb.TB),1);
        for i = 1:length(perturb.TB)
            muTB(i) = cspice_bodvrd(perturb.TB{i},'GM',1);
        end

        switch ref_sys.obs
            
            case 'EARTH'
                % check gravity degree
                if perturb.n > n_hifi(1) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

                filename = 'arael/data/hifi/Earth/EGM96.txt';

            case 'MOON'
                % check gravity degree
                if perturb.n > n_hifi(2) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

                filename = 'arael/data/hifi/Moon/gggrx_1200a_sha.txt';

            case 'MARS'
                % check gravity degree
                if perturb.n > n_hifi(3) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

                filename = 'arael/data/hifi/Mars/GGM1025A.txt';

            case 'VENUS'
                % check gravity degree
                if perturb.n > n_hifi(4) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

                filename = 'arael/data/hifi/Venus/shgj180u.txt';

        end

        % compute expansion coefficients
        fprintf('Computing expansion coefficients...\n')
        [Cnm_mod, Snm_mod] = gravCoeff(filename,perturb.n);

        % Compute legendre polynomial and derivative as functions
        fprintf('Computing harmonic functions...\n');
        Pnm_norm = expansionFunc(perturb.n);
      
        % GRAVITY
        g = @(t,r) gravAcc(t,r,mu,Rp,Pnm_norm,Cnm_mod,Snm_mod,perturb.n,ref_sys.inertial,ref_sys.obs,et);

        % THIRD-BODY
        aTB = @(t,r) pertTB(t,r,perturb.TB,muTB,et,ref_sys.inertial,ref_sys.obs);

        % SOLAR RADIATION PRESSURE
        aSRP = @(t,r) pertSRP(t,r,Rp,perturb.SRP,spacecraft.cR,spacecraft.A,spacecraft.m,ref_sys.inertial,ref_sys.obs,et);

        % PERTURBING ACCELERATION
        aTOT = @(t,x) g(t,x(1:3)) + aTB(t,x(1:3)) + aSRP(t,x(1:3));

        % INITIAL CONDITION
        equi0 = car2equi(init_cond.x0,mu);

        % integrate using Equinoctial elements
        fprintf('Propagating the orbit...\n');
        tic
        [t,equi_s] = ode113(@(t,equi) hifi_rhs(t,equi,aTOT,mu),init_cond.tSpan,equi0,options);
        toc

        % convert equinoctial elements into state
        y = zeros(length(t),6);
        for i = 1:length(t)
            y(i,:) = equi2car(equi_s(i,:),mu);
        end


    case 'approx'
        % check for allowed observers
        obs_allowed_approx = {'EARTH';'MOON'};

        if ismember(ref_sys.obs,obs_allowed_approx) ~= 1
            fprintf('\n')
            fprintf('ERROR: invalid observer body \n')
            return
        end

        switch ref_sys.obs

            case 'EARTH'
                % check gravity degree
                if perturb.n > n_approx(1) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

            case 'MOON'
                % check gravity degree
                if perturb.n > n_approx(2) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

            case 'MARS'
                % check gravity degree
                if perturb.n > n_approx(3) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end


            case 'VENUS'
                % check gravity degree
                if perturb.n > n_approx(4) || perturb.n < 0
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

        end

        % retrieve the zonal harmonics coeffieicnts
        switch(ref_sys.obs)

            case 'EARTH'
                JN = [0.00108263,-2.5321530e-6,-1.6109877e-6,-2.3578565e-7,5.4316985e-7];

            case 'MOON'
                JN = [2.0323e-4,8.4759e-06,-9.5919e-06];

            case 'MARS'
                JN = [1.9555e-03,3.1450e-05,-1.5377e-05];

            case 'VENUS'
                JN = [4.4044e-06,-2.1082e-06,-2.1474e-06];

        end

        % compute the gravitational constant and radius
        mu = cspice_bodvrd(ref_sys.obs,'GM',1);
        Rp = mean(cspice_bodvrd(ref_sys.obs,'RADII',3));
        
        muTB = zeros(length(perturb.TB),1);
        for i = 1:length(perturb.TB)
            muTB(i) = cspice_bodvrd(perturb.TB{i},'GM',1);
        end

        % GRAVITY
        aZH = @(t,r) pertZH(t,r,ref_sys.inertial,ref_sys.obs,mu,Rp,perturb.n,JN,et);

        % THIRD-BODY
        aTB = @(t,r) pertTB(t,r,perturb.TB,muTB,et,ref_sys.inertial,ref_sys.obs);

        % SOLAR RADIATION PRESSURE
        aSRP = @(t,r) pertSRP(t,r,Rp,perturb.SRP,spacecraft.cR,spacecraft.A,spacecraft.m,ref_sys.inertial,ref_sys.obs,et);

        % PERTURBING ACCELERATION
        aTOT = @(t,x) aZH(t,x(1:3)) + aTB(t,x(1:3)) + aSRP(t,x(1:3));

        % INITIAL CONDITION
        equi0 = car2equi(init_cond.x0,mu);

        % integrate using Equinoctial elements
        fprintf('Propagating the orbit...\n');
        tic
        [t,equi_s] = ode113(@(t,equi) approx_rhs(t,equi,aTOT,mu),init_cond.tSpan,equi0,options);
        toc

        % convert equinoctial elements into state
        y = zeros(length(t),6);
        for i = 1:length(t)
            y(i,:) = equi2car(equi_s(i,:),mu);
        end

    case 'averaged'

    case 'full'
        % check for allowed observers
        obs_allowed_full = {'SUN'};

        if ismember(ref_sys.obs,obs_allowed_full) ~= 1
            fprintf('\n')
            fprintf('ERROR: invalid observer body \n')
            return
        end

        % compute the gravitational constantants
        mu = cspice_bodvrd(ref_sys.obs,'GM',1);
        
        muTB = zeros(length(perturb.TB),1);
        for i = 1:zeros(length(perturb.TB),1)
            muTB(i) = cspice_bodvrd(perturb.TB{i},'GM',1);
        end

        % SOLAR RADIATION PRESSURE
        aSRP = @(t,r) pertSRP_full(t,r,perturb.SRP,spacecraft.cR,spacecraft.A,spacecraft.m,ref_sys.inertial,ref_sys.obs,et);

        % PERTURBING ACCELERATION
        aTOT = @(t,x) aSRP(t,x(1:3));

        % integrate using cartesian state
        fprintf('Propagating the orbit...\n');
        tic
        [t,y] = ode113(@(t,x) full_rhs(t,x,mu,perturb.TB,muTB,aTOT,ref_sys.inertial,ref_sys.obs,et),init_cond.tSpan,init_cond.x0,options);
        toc
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ODE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------- HIFI ------------------------------------%

function dequi = hifi_rhs(t,equi,acc_ijk,mu)

% convert state into cartesian elements
x = equi2car(equi,mu);

% convert acceleration in rtn
r_ijk = x(1:3);
v_ijk = x(4:6);
r_vers = r_ijk/norm(r_ijk);
rv_cross = cross(r_ijk, v_ijk);
n_vers = rv_cross/norm(rv_cross);
t_vers = cross(n_vers, r_vers);
IJK_2_RTN = [r_vers, t_vers, n_vers]';
acc_rtn = IJK_2_RTN * acc_ijk(t,x);

% Compute derivative in equinoctial elements
dequi = gaussEquinoctial(equi, mu, acc_rtn);

end

%------------------------------ APPROX -----------------------------------%

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
IJK_2_RTN = [r_vers, t_vers, n_vers]';
acc_rtn = IJK_2_RTN * acc_ijk(t,x);

% compute variation of equinoctial elements
dequi = gaussEquinoctial(equi,mu,acc_rtn);

end

%------------------------------- FULL ------------------------------------%

function [dxdt] = full_rhs(t,x,mu,NB,muNB,acc_ijk,ref_sys,obs,et)

% Initialize right-hand-side
dxdt = zeros(6,1);

% Position detivative is object's velocity
dxdt(1:3) = x(4:6);

% Extract the object position from state x
rSC = x(1:3);

% Compute acceleration due to main attractor
dxdt(4:6) = -mu*rSC/(norm(rSC)^3);

% Loop over all bodies
for i = 1:size(NB(:,1))

    % compute n-bodies position and gravitational constant
    rNB = cspice_spkpos(NB{i},t+et,ref_sys,'NONE',obs);

    % Extract object position wrt. i-th celestial body
    rRel = rSC - rNB;

    % Compute square distance and distance
    rRelNorm = norm(rRel);

    % Compute the gravitational acceleration using Newton's law
    acc_pos_i = -muNB(i)*rRel/rRelNorm^3;

    % Sum up acceleration to right-hand-side
    dxdt(4:6) = dxdt(4:6) + acc_pos_i;
end

dxdt(4:6) = dxdt(4:6) + acc_ijk(t,x);

end