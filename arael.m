function [t,y] = arael(init_cond,ref_sys,perturb,settings)
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
%  - SOLAR RADIATION PRESSURE: not yet implemented
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
%  - SOLAR RADIATION PRESSURE: not yet implemented
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
%  - SOLAR RADIATION PRESSURE: not yet implemented
%
%-------------------------------------------------------------------------%
% THIS FUNCTION MAKE USE OF THE SPICE TOOLKIT DEVELOPED BY NASA, LEARN
% MORE ABOUT SPICE AT: https://naif.jpl.nasa.gov/naif/toolkit.html
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
% TO DO: - debug MARS and VENUS using GMAT
%        - check zonal coefficients for MOON,MARS,VENUS
%        - add 'average'
%        - add SRP and DRAG
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
obs_allowed = {'SUN';'EARTH';'MOON';'MARS';'VENUS'};

% gravity
n_hifi = [360;    % Earth
    1200;         % Moon
    80;           % Mars
    180];         % Venus      

n_approx = [6;    % Earth
       2];        % Moon

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
                if perturb.n > n_hifi(1)
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

                filename = 'arael/data/hifi/Earth/EMG96.txt';

            case 'MOON'
                % check gravity degree
                if perturb.n > n_hifi(2)
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

                filename = 'arael/data/hifi/Moon/gggrx_1200a_sha.txt';

            case 'MARS'
                % check gravity degree
                if perturb.n > n_hifi(3)
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

                filename = 'arael/data/hifi/Mars/GGM1025A.txt';

            case 'VENUS'
                % check gravity degree
                if perturb.n > n_hifi(4)
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
        g = @(t,r) gravAcc(t,r,mu,Rp,Pnm_norm,Cnm_mod,Snm_mod,perturb.n,ref_sys.inertial,ref_sys.obs,init_cond.et);

        % THIRD-BODY
        aTB = @(t,r) pertTB(t,r,perturb.TB,muTB,init_cond.et,ref_sys.inertial,ref_sys.obs);

        % PERTURBING ACCELERATION
        aTOT = @(t,x) g(t,x(1:3)) + aTB(t,x(1:3));

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
                if perturb.n > n_approx(1)
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end

            case 'MOON'
                % check gravity degree
                if perturb.n > n_approx(2)
                    fprintf('ERROR: invalid order of the gravitational field')
                    return
                end
        end

        % retrieve the zonal harmonics coeffieicnts
        switch(obs)

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
        aZH = @(t,r) pertZH(t,r,ref_sys.inertial,ref_sys.obs,mu,Rp,perturb.n,JN,init_cond.et);

        % THIRD-BODY
        aTB = @(t,r) pertTB(t,r,perturb.TB,muTB,init_cond.et,ref_sys.inertial,ref_sys.obs);

        % PERTURBING ACCELERATION
        aTOT = @(t,x) aZH(t,x(1:3)) + aTB(t,x(1:3));

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

        % integrate using cartesian state
        fprintf('Propagating the orbit...\n');
        tic
        [t,y] = ode113(@(t,x) full_rhs(t,x,mu,perturb.TB,muTB,ref_sys.inertial,ref_sys.obs,init_cond.et),init_cond.tSpan,init_cond.x0,options);
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
IJK_2_RTN = [r_vers, t_vers, n_vers];
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
IJK_2_RTN = [r_vers, t_vers, n_vers];
acc_rtn = IJK_2_RTN * acc_ijk(t,x);

% compute variation of equinoctial elements
dequi = gaussEquinoctial(equi,mu,acc_rtn);

end

%------------------------------- FULL ------------------------------------%

function [dxdt] = full_rhs(t,x,mu,NB,muNB,ref_sys,obs,et)

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

end