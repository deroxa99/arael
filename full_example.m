%%% ARAEL TEST
clearvars
close all
clc

%% Data

init_cond.x0 = cspice_spkezr('EARTH',4e06,'ECLIPJ2000','NONE','SUN');

% time
init_cond.et = 0;
init_cond.tSpan = 0:60:365*24*3600;

% gravity
perturb.n = 0;

% third body
%perturb.TB = {'EARTH'};
perturb.TB = {'MERCURY';
    'VENUS';
    'EARTH';
    'MOON';
    'MARS BARYCENTER';
    'JUPITER BARYCENTER';
    'SATURN BARYCENTER';
    'URANUS BARYCENTER';
    'NEPTUNE BARYCENTER';
    'PLUTO BARYCENTER'};

% ref sys
ref_sys.inertial = 'ECLIPJ2000';
%ref_sys.inertial = 'J2000';
ref_sys.obs = 'SUN';

% settings
settings.mode = 'full';
settings.rel_tol = 1e-09;
settings.abs_tol = 1e-10;

%% integrate
[t,y] = arael(init_cond,ref_sys,perturb,settings);

%% post processing

% 3D plot
figure(1)
plot3(y(:,1),y(:,2),y(:,3))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal