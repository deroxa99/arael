%%% ARAEL TEST
clearvars
close all
clc

%% Data

% initial condition
a = 4.214988153818061;
e = 0.000002586353876;
i = 0.000022934070698;
OM = 0.000022941161937;
om = 0.000550639719503;
f = 0.000414587334010;
kep0 = [a,e,i,OM,om,f]*1e04;

%mu = 4.9028e+03; % Moon
mu = 3.9860e+05; % Earth

init_cond.x0 = kep2car(kep0,mu);

% time
init_cond.utc0 ='2024-02-05 05:17:42.819936 UTC';
init_cond.tSpan = 0:60:10*24*3600;

% gravity
perturb.n = 10;

% third body
%perturb.TB = {'EARTH'};
perturb.TB = {'SUN','MOON'};

% ref sys
ref_sys.inertial = 'J2000';
%ref_sys.inertial = 'J2000';
ref_sys.obs = 'EARTH';

% settings
settings.mode = 'hifi';
settings.rel_tol = 1e-09;
settings.abs_tol = 1e-10;

% satellite

%% integrate
[t,y] = arael(init_cond,ref_sys,perturb,settings);

%% post processing

% retreieve keplerians
kep = zeros(length(t),6);

for i = 1:length(t)
    kep(i,:) = car2kep(y(i,:),mu);
end

%% plot

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
earth_J6 = 'arael-main\utils\gmat_debug\EARTH_debug.txt';
earth_3B_J6 = 'arael-main\utils\gmat_debug\EARTH_debug_3B_J6.txt';
earth_3B = 'arael-main\utils\gmat_debug\EARTH_debug_3B.txt';
moon_J2 = 'arael-main\utils\gmat_debug\MOON_debug.txt';
moon_3B = 'arael-main\utils\gmat_debug\MOON_debug_3B.txt';
moon_3B_E = 'arael-main\utils\gmat_debug\MOON_debug_3B_E.txt';
earth_sgp4 = 'arael-main\utils\gmat_debug\EARTH_debug_sgp4.txt';

data = readmatrix(earth_sgp4);
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
plot3(y_gmat(:,1),y_gmat(:,2),y_gmat(:,3),'b')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')


% plot keplerians
figure(4)
sgtitle('Evolution of keplerian elements')
subplot(2,3,1)
hold on
plot(t, kep(:,1),'r',t_gmat,kep_gmat(:,1),'b')
title('Semi-major axis')
xlabel('t [Orbits]')
ylabel('a [km]')
grid on
subplot(2,3,2)
hold on
plot(t, kep(:,2),'r',t_gmat,kep_gmat(:,2),'b')
title('Eccectricity')
xlabel('t [Orbits]')
ylabel('e []')
grid on
subplot(2,3,3)
hold on
plot(t, kep(:,3),'r',t_gmat,kep_gmat(:,3),'b')
title('Inclination')
xlabel('t [Orbits]')
ylabel('i [°]')
grid on
subplot(2,3,4)
hold on
plot(t, kep(:,4),'r',t_gmat,kep_gmat(:,4),'b')
title('Right Ascension of the ascending node')
xlabel('t [Orbits]')
ylabel('\Omega [°]')
grid on
subplot(2,3,5)
hold on
plot(t, kep(:,5),'r',t_gmat,kep_gmat(:,5),'b')
title('Argument of Pericenter')
xlabel('t [Orbits]')
ylabel('\omega [°]')
grid on
subplot(2,3,6)
hold on
plot(t, wrapTo360(kep(:,6)),'r',t_gmat,kep_gmat(:,6),'b')
title('True Anomaly')
xlabel('t [Orbits]')
ylabel('\theta [°]')
grid on

%% compare with sgp4

%%% sgp4
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)
et0 = cspice_str2et(init_cond.utc0);

% load the tle (breeze 26373)
tle_1 = '1 26373U 00029B   24036.22063449 -.00000084  00000-0  00000-0 0  9990';
tle_2 = '2 26373  13.1105  12.8500 0258164 316.2858 239.8887  1.00328814 86779';

% compute state
satrec = twoline2rv(tle_1, tle_2, typerun,'e', opsmode, whichconst);

% retrieve tle epoch
[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% integrate using sgp4
tSpan_sgp4 = (init_cond.tSpan + et0 - sat_epoch_et)/60;
r_teme = zeros(3, length(tSpan_sgp4));
v_teme = zeros(3, length(tSpan_sgp4));

for i = 1:length(tSpan_sgp4) 
[~,r_teme(:,i),v_teme(:,i)] = sgp4(satrec,  tSpan_sgp4(i));
end

% convert from TEME to ECI
ateme = [0;0;0];
ttt = cspice_unitim(init_cond.tSpan + et0, 'ET', 'TDT')/cspice_jyear()/100;
arcsec2rad = pi/(180*3600);
dPsi = -0.113638*arcsec2rad;
dEps = -0.007048*arcsec2rad;
r_sgp4 = zeros(3, length(tSpan_sgp4));
v_sgp4 = zeros(3, length(tSpan_sgp4));

for i = 1:length(tSpan_sgp4) 
[r_sgp4(:,i), v_sgp4(:,i), ~] = teme2eci(r_teme(:,i), v_teme(:,i), ateme, ttt(i), dPsi, dEps);
end

y_sgp4 = [r_sgp4;v_sgp4]';
t_sgp4 = init_cond.tSpan;

%% plot

% retreieve keplerians
kep_sgp4 = zeros(length(t),6);

for i = 1:length(t)
    kep_sgp4(i,:) = car2kep(y_sgp4(i,:),mu);
end

% 3D plot
figure(3)
plot3(y(:,1),y(:,2),y(:,3),'r')
hold on
plot3(y_sgp4(:,1),y_sgp4(:,2),y_sgp4(:,3),'g')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')


% plot keplerians
figure(4)
sgtitle('Evolution of keplerian elements')
subplot(2,3,1)
hold on
plot(t, kep(:,1),'r',t_sgp4,kep_sgp4(:,1),'g')
title('Semi-major axis')
xlabel('t [Orbits]')
ylabel('a [km]')
grid on
subplot(2,3,2)
hold on
plot(t, kep(:,2),'r',t_sgp4,kep_sgp4(:,2),'g')
title('Eccectricity')
xlabel('t [Orbits]')
ylabel('e []')
grid on
subplot(2,3,3)
hold on
plot(t, kep(:,3),'r',t_sgp4,kep_sgp4(:,3),'g')
title('Inclination')
xlabel('t [Orbits]')
ylabel('i [°]')
grid on
subplot(2,3,4)
hold on
plot(t, kep(:,4),'r',t_sgp4,kep_sgp4(:,4),'g')
title('Right Ascension of the ascending node')
xlabel('t [Orbits]')
ylabel('\Omega [°]')
grid on
subplot(2,3,5)
hold on
plot(t, kep(:,5),'r',t_sgp4,kep_sgp4(:,5),'g')
title('Argument of Pericenter')
xlabel('t [Orbits]')
ylabel('\omega [°]')
grid on
subplot(2,3,6)
hold on
plot(t, wrapTo360(kep(:,6)),'r',t_sgp4,kep_sgp4(:,6),'g')
title('True Anomaly')
xlabel('t [Orbits]')
ylabel('\theta [°]')
grid on
