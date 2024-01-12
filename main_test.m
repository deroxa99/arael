%%% ARAEL TEST
clearvars
close all
clc

%% Data

% initial condition
a = 7000;
e = 0;
i = 0;
OM = 0;
om = 0;
f = 0;
kep0 = [a,e,i,OM,om,f];

%mu = 4.9028e+03; % Moon
mu = 3.9860e+05; % Earth

init_cond.x0 = kep2car(kep0,mu);

% time
init_cond.et = 0;
init_cond.tSpan = 0:60:10*24*3600;

% gravity
perturb.n = 0;

% third body
%perturb.TB = {'EARTH'};
perturb.TB = {'SUN','MOON'};

% ref sys
ref_sys.inertial = 'J2000';
%ref_sys.inertial = 'J2000';
ref_sys.obs = 'EARTH';

% settings
settings.mode = 'approx';
settings.rel_tol = 1e-09;
settings.abs_tol = 1e-10;

%% integrate
[t,y] = arael(init_cond,ref_sys,perturb,settings);

%% post processing

% retreieve keplerians
kep = zeros(length(t),6);

for i = 1:length(t)
    kep(i,:) = car2kep(y(i,:),mu);
end

%%% plot

% 3D plot
figure(1)
plot3(y(:,1),y(:,2),y(:,3))
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')


figure(2)
sgtitle('Evolution of keplerian elements')
subplot(2,3,1)
hold on
plot(t, kep(:,1),'r')
title('Semi-major axis')
xlabel('t [Orbits]')
ylabel('a [km]')
grid on
subplot(2,3,2)
hold on
plot(t, kep(:,2),'r')
title('Eccectricity')
xlabel('t [Orbits]')
ylabel('e []')
grid on
subplot(2,3,3)
hold on
plot(t, kep(:,3),'r')
title('Inclination')
xlabel('t [Orbits]')
ylabel('i [°]')
grid on
subplot(2,3,4)
hold on
plot(t, kep(:,4),'r')
title('Right Ascension of the ascending node')
xlabel('t [Orbits]')
ylabel('\Omega [°]')
grid on
subplot(2,3,5)
hold on
plot(t, kep(:,5),'r')
title('Argument of Pericenter')
xlabel('t [Orbits]')
ylabel('\omega [°]')
grid on
subplot(2,3,6)
hold on
plot(t, wrapTo360(kep(:,6)),'r')
title('True Anomaly')
xlabel('t [Orbits]')
ylabel('\theta [°]')
grid on

%% comparison with gmat
earth_J6 = 'arael\utils\gmat_debug\EARTH_debug.txt';
earth_3B_J6 = 'arael\utils\gmat_debug\EARTH_debug_3B_J6.txt';
earth_3B = 'arael\utils\gmat_debug\EARTH_debug_3B.txt';
moon_J2 = 'arael\utils\gmat_debug\MOON_debug.txt';
moon_3B = 'arael\utils\gmat_debug\MOON_debug_3B.txt';
moon_3B_E = 'arael\utils\gmat_debug\MOON_debug_3B_E.txt';

data = readmatrix(earth_3B);
y_gmat = data(:,2:7);
t_gmat = data(:,8);

err = vecnorm(y - y_gmat)

% retreieve keplerians
kep_gmat = zeros(length(t),6);

for i = 1:length(t)
    kep_gmat(i,:) = car2kep(y_gmat(i,:),mu);
end

% 3D plot
figure(3)
plot3(y(:,1),y(:,2),y(:,3),'r')
hold on
plot3(y_gmat(:,1),y_gmat(:,2),y_gmat(:,3),'g')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')


% plot keplerians
figure(4)
sgtitle('Evolution of keplerian elements')
subplot(2,3,1)
hold on
plot(t, kep(:,1),'r',t_gmat,kep_gmat(:,1),'g')
title('Semi-major axis')
xlabel('t [Orbits]')
ylabel('a [km]')
grid on
subplot(2,3,2)
hold on
plot(t, kep(:,2),'r',t_gmat,kep_gmat(:,2),'g')
title('Eccectricity')
xlabel('t [Orbits]')
ylabel('e []')
grid on
subplot(2,3,3)
hold on
plot(t, kep(:,3),'r',t_gmat,kep_gmat(:,3),'g')
title('Inclination')
xlabel('t [Orbits]')
ylabel('i [°]')
grid on
subplot(2,3,4)
hold on
plot(t, kep(:,4),'r',t_gmat,kep_gmat(:,4),'g')
title('Right Ascension of the ascending node')
xlabel('t [Orbits]')
ylabel('\Omega [°]')
grid on
subplot(2,3,5)
hold on
plot(t, kep(:,5),'r',t_gmat,kep_gmat(:,5),'g')
title('Argument of Pericenter')
xlabel('t [Orbits]')
ylabel('\omega [°]')
grid on
subplot(2,3,6)
hold on
plot(t, wrapTo360(kep(:,6)),'r',t_gmat,kep_gmat(:,6),'g')
title('True Anomaly')
xlabel('t [Orbits]')
ylabel('\theta [°]')
grid on