clearvars; close all; clc;


% SPHERICAL MOTION WITH GG only PERTURBATIONS in GEO. 
% FREE DRIFT + CONTROLLED MOTION 


%% SETUP
hours = 3600;                                    %Hours to seconds
days = 24*hours;                                 %Days to seconds


ll = 1e5; tt = 1e5;                            %Uncomment to adimensionalize
dat.Re = 6378.137; %/ll; 
dat.mue = 398600.4415; %*tt^2/ll^3;
dat.Rgeo = 42165.8; %/ll;                             
dat.w = sqrt(dat.mue/(dat.Rgeo)^3);

% values taken from 'vallado' 
grav.c20 = -1.083e-3;  
grav.c21 = -2.414e-10;  
grav.c22 = 1.574e-6;  
grav.c30 = 2.532e-6;  
grav.c31 = 2.191e-6; 
grav.c32 = 3.089e-7;  
grav.c33 = 1.006e-7;  
grav.s21 = 1.543e-9;  
grav.s22 = -9.038e-7;  
grav.s31 = 2.687e-7; 
grav.s32 = -2.115e-7;
grav.s33 = 1.972e-7;

% Set initial FREE-DRIFT conditions and times
l_nom = 60;                           %nominal longitude #degrees
r0 = dat.Rgeo;                        %initial radius (altitude = Rgeo - Re)
l0 = deg2rad(l_nom - 0.04);           %initial longitude error 
phi0 = 0;                             %initial latitude 

y0_FD = [r0, l0, phi0, 0, 0, 0]; 
% writematrix(y0_FD, 'init.txt');


% t0 = date2mjd2000([2000, 12, 21, 0, 0, 0])*days;                           % set initial date and convert it in seconds
t0 = 0; 
tf1 = t0 + 2.5*days;                                                       % final free drift time (in seconds) 
tf2 = tf1 + 0.5*days;                                                      % final control time (in seconds)
nit = 1000; 
t_FREEDRIFT = linspace(t0, tf1, nit); %./tt;
t_CONTROL = linspace(tf1, tf2, nit);  %./tt; 

%% SPHERICAL MOTION -- FREE DRIFT

% long = atan2(y, x);
% lat = atan2(z, sqrt(x.^2 + y.^2));
% r = sqrt(x.^2 + y.^2 + z.^2);

% Use ODE113 to integrate spherical equations from t0 to tf:

options = odeset('reltol', 1.e-12,'abstol', 1.e-12);
[t_FD, sph_FD] = ode78(@(t,sph) SPHmotionGG_FREEDRIFT(t, sph, dat, grav), t_FREEDRIFT, y0_FD, options);

% Save txt file for python DA initial condition 
final_FD = sph_FD(end,:)';
writematrix(final_FD, 'fin_FREEDRIFT.txt'); 

%% DIMENSIONALIZE 
% final_FD(1) = final_FD(1) * ll;
% final_FD(5) = final_FD(5) / tt;
% final_FD(6) = final_FD(6) / tt;
% 
% dat.Re = 6378.137; 
% dat.mue = 398600.4415;
% dat.Rgeo = 42165.8;                              %nominal radius of geostationary orbit
% dat.w = sqrt(dat.mue/(dat.Rgeo)^3);


%% COMPUTE SPHERICAL COSTATE DYNAMICS WITH SYMBOLIC FOR CONTROL
% [cost_dyn] = costate; 

%% SPHERICAL MOTION -- CONTROL -- SETUP conditions  

% Set initial control conditions----State and Costate 
y0_C_ref = zeros(6,1);     %Initial control state is final free drift state
deltaY0_C = final_FD;                 
y0_C = y0_C_ref + deltaY0_C; 

lamda0_ref = zeros(6,1); 
Dlamda0 = [-2e-12; -1e-12; 0; 2e-7; 1e-7; 0];
% Dlamda0 = importdata('LAMBDA0.txt');
lamda0 = lamda0_ref + Dlamda0; 


dyn0 = [y0_C; lamda0];

%% SPHERICAL MOTION -- CONTROL -- PROPAGATIONS. Spherical equations from t0 to tf:

options = odeset('reltol', 1.e-12,'abstol', 1.e-12);
[t_C, sph_C] = ode78(@(t,sph) SPHmotionGG_CONTROL(t, sph, dat, grav), t_CONTROL, dyn0, options);

writematrix(sph_C(end,:)', 'fin_CONTROL.txt');

%% SPHERICAL MOTION -- CONTROL -- SHOOTING

state.initial = final_FD;            %end of free drift
state.final = y0_FD';                %Back to initial free drift state for SK

% fsolve options
options_fsolve = optimoptions(@fsolve, 'Display', 'iter','Maxiter', 1000,'FunctionTolerance', 1e-12, 'StepTolerance',1e-6);
% fsolve call-- O = output, contains final time and initial conditions for
% lambda_r, lambda_v
O = fsolve(@opt_solv, lamda0, options_fsolve, state, t_CONTROL, dat, grav);

dyn0_shoot = [state.initial; O]; 

% Propagation to compute state and costate evolution
[t_S, sph_S] = ode78(@(t,sph)  SPHmotionGG_CONTROL(t, sph, dat, grav), t_CONTROL, dyn0_shoot, options);


%% PLOT SPHERICAL free drift + control

figure()

subplot(3,1,1)
plot((t_FD - t0)/24/3600, sph_FD(:, 1))
title('Distance')
xlabel('Days')
grid on
grid minor
axis tight

subplot(3,1,2)
plot((t_FD - t0)/24/3600, ((rad2deg(sph_FD(:, 2)) - l_nom)))
title('Longitude error')
xlabel('Days')
grid on
grid minor
axis tight

subplot(3,1,3)
plot((t_FD - t0)/24/3600, rad2deg(sph_FD(:, 3)))
title('Latitude')
xlabel('Days')
grid on
grid minor
axis tight


figure()  % Groundtrack
plot(rad2deg(sph_FD(:,2)), rad2deg(sph_FD(:,3)), '-r','LineWidth',2)
hold on
grid on
title('Groundtrack');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
plot(rad2deg(sph_FD(1,2)), rad2deg(sph_FD(1,3)),'sm','LineWidth',1.5); 
plot(rad2deg(sph_FD(end,2)), rad2deg(sph_FD(end,3)),'ob','LineWidth',2);

plot(rad2deg(sph_C(:,2)), rad2deg(sph_C(:,3)), '--g','LineWidth',1.5)
plot(rad2deg(sph_C(1,2)), rad2deg(sph_C(1,3)),'sg','LineWidth', 2); 
plot(rad2deg(sph_C(end,2)), rad2deg(sph_C(end,3)),'ok','LineWidth',2);
legend('Track', 'Start Free drift', 'End Free drift', 'Controlled track', 'Start control', 'End control', Location='best')

%% SHOOTING PLOT

figure()  % Groundtrack
plot(rad2deg(sph_FD(:,2)), rad2deg(sph_FD(:,3)), '-.k','LineWidth',2)
hold on
grid on
title('Groundtrack');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
plot(rad2deg(sph_FD(1,2)), rad2deg(sph_FD(1,3)),'sm','LineWidth',1.5); 
plot(rad2deg(sph_FD(end,2)), rad2deg(sph_FD(end,3)),'ob','LineWidth',2);

plot(rad2deg(sph_S(:,2)), rad2deg(sph_S(:,3)), '-g','LineWidth',1.5)
plot(rad2deg(sph_S(1,2)), rad2deg(sph_S(1,3)),'sg','LineWidth', 2); 
plot(rad2deg(sph_S(end,2)), rad2deg(sph_S(end,3)),'ok','LineWidth',2);
legend('Track', 'Start Free drift', 'End Free drift', 'Controlled track', 'Start control', 'End control', Location='best')

% Control evolution

ur = - sph_S(:,10); 

for i = 1 : length(t_S)
    ul(i) = - sph_S(i,11) / (sph_S(i,1) * cos(sph_S(i,3)));
    uphi(i) = - sph_S(i,12) / sph_S(i,1);
end

figure()
plot((t_S - t0)/24/3600, ur, 'k')
hold on
plot((t_S - t0)/24/3600, ul, '--k')
plot((t_S - t0)/24/3600, uphi, '-.k')
legend('u_{r}', 'u_{l}', 'u_{\phi}', Location='best')
xlabel('Days')
ylabel('m/s^{2}')

%% FUNCTIONS 

function dsph = SPHmotionGG_FREEDRIFT(~, sph, dat, grav)

    mu = dat.mue; omega = dat.w;  R = dat.Re; 
    

    r = sph(1);             v = sph(4); 
    l = sph(2);             csi = sph(5);
    phi = sph(3);           n = sph(6);
    
    c20 = grav.c20;  c21 = grav.c21;  c22 = grav.c22;  c30 = grav.c30;  c31 = grav.c31; 
    c32 = grav.c32;  c33 = grav.c33;  
    s21 = grav.s21;  s22 = grav.s22;  s31 = grav.s31;  s32 = grav.s32;  s33 = grav.s33;


    % GRAVITY PERTURBING ACCELERATIONS (zonal + tesseral) --- derivatives taken with symbolic 

    aPg_r = -mu*((3*R^2*c20*(3*sin(phi)^2 - 1))/(2*r^4) + ...
            (9*R^2*cos(phi)^2*(s22*sin(2*l) + c22*cos(2*l)))/r^4 + ...
            (60*R^3*cos(phi)^3*(s33*sin(3*l) + c33*cos(3*l)))/r^5 + ...
            (2*R^3*c30*sin(phi)*(5*sin(phi)^2 - 3))/r^5 + ...
            (9*R^2*cos(phi)*sin(phi)*(c21*cos(l) + s21*sin(l)))/r^4 + ...
            (60*R^3*cos(phi)^2*sin(phi)*(s32*sin(2*l) + c32*cos(2*l)))/r^5 + ...
            (2*R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/r^5);
 


 
    aPg_l =  -mu*((3*R^2*cos(phi)^2*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^3 + ...
             (15*R^3*cos(phi)^3*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^4 + ...
             (3*R^2*cos(phi)*sin(phi)*(c21*sin(l) - s21*cos(l)))/r^3 + ...
             (15*R^3*cos(phi)^2*sin(phi)*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + ...
             (R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/(2*r^4));
 

 
    aPg_phi = mu*((15*R^3*cos(phi)^3*(s32*sin(2*l) + c32*cos(2*l)))/r^4 + ...
              (3*R^2*cos(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 - ...
              (3*R^2*sin(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 + ...
              (R^3*c30*cos(phi)*(5*sin(phi)^2 - 3))/(2*r^4) - ...
              (30*R^3*cos(phi)*sin(phi)^2*(s32*sin(2*l) + c32*cos(2*l)))/r^4 - ...
              (45*R^3*cos(phi)^2*sin(phi)*(s33*sin(3*l) + c33*cos(3*l)))/r^4 + ...
              (5*R^3*c30*cos(phi)*sin(phi)^2)/r^4 - ...
              (R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/(2*r^4) + ...
              (15*R^3*cos(phi)^2*sin(phi)*(c31*cos(l) + s31*sin(l)))/r^4 - ...
              (6*R^2*cos(phi)*sin(phi)*(s22*sin(2*l) + c22*cos(2*l)))/r^3 + ...
              (3*R^2*c20*cos(phi)*sin(phi))/r^3);
 

    % FINAL DYNAMICS
    dr = v; 
    dl = csi; 
    dphi = n; 
    dv = -mu/r^2 + r*n^2 + r*(csi + omega)^2 *(cos(phi))^2 + ...
          aPg_r; 

    dcsi = 2*n*(csi + omega)*tan(phi) - 2*v/r *(csi + omega) + ...
           aPg_l/(r * cos(phi))^2; 

    dn = -2*v/r *n - (csi + omega)^2 *sin(phi)*cos(phi) + ...
          aPg_phi/r^2; 


    dsph = [dr; dl; dphi; dv; dcsi; dn]; 
end



function dsph = SPHmotionGG_CONTROL(~, dyn, dat, grav)

    mue = dat.mue; om = dat.w;  R = dat.Re;     

    r = dyn(1);             v = dyn(4); 
    l = dyn(2);             csi = dyn(5);
    phi = dyn(3);           n = dyn(6);
    l_r = dyn(7);           l_v = dyn(10);
    l_l = dyn(8);           l_csi = dyn(11);
    l_phi = dyn(9);         l_n = dyn(12);
                                  
    c20 = grav.c20;  c21 = grav.c21;  c22 = grav.c22;  c30 = grav.c30;  c31 = grav.c31; 
    c32 = grav.c32;  c33 = grav.c33;  
    s21 = grav.s21;  s22 = grav.s22;  s31 = grav.s31;  s32 = grav.s32;  s33 = grav.s33;

 
    % GRAVITY PERTURBING ACCELERATIONS (zonal + tesseral) --- derivatives taken with symbolic 
    aPg_r = -mue*((3*R^2*c20*(3*sin(phi)^2 - 1))/(2*r^4) + ...
            (9*R^2*cos(phi)^2*(s22*sin(2*l) + c22*cos(2*l)))/r^4 + ...
            (60*R^3*cos(phi)^3*(s33*sin(3*l) + c33*cos(3*l)))/r^5 + ...
            (2*R^3*c30*sin(phi)*(5*sin(phi)^2 - 3))/r^5 + ...
            (9*R^2*cos(phi)*sin(phi)*(c21*cos(l) + s21*sin(l)))/r^4 + ...
            (60*R^3*cos(phi)^2*sin(phi)*(s32*sin(2*l) + c32*cos(2*l)))/r^5 + ...
            (2*R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/r^5);
 

    aPg_l =  -mue*((3*R^2*cos(phi)^2*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^3 + ...
             (15*R^3*cos(phi)^3*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^4 + ...
             (3*R^2*cos(phi)*sin(phi)*(c21*sin(l) - s21*cos(l)))/r^3 + ...
             (15*R^3*cos(phi)^2*sin(phi)*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + ...
             (R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/(2*r^4));
 
 
    aPg_phi = mue*((15*R^3*cos(phi)^3*(s32*sin(2*l) + c32*cos(2*l)))/r^4 + ...
              (3*R^2*cos(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 - ...
              (3*R^2*sin(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 + ...
              (R^3*c30*cos(phi)*(5*sin(phi)^2 - 3))/(2*r^4) - ...
              (30*R^3*cos(phi)*sin(phi)^2*(s32*sin(2*l) + c32*cos(2*l)))/r^4 - ...
              (45*R^3*cos(phi)^2*sin(phi)*(s33*sin(3*l) + c33*cos(3*l)))/r^4 + ...
              (5*R^3*c30*cos(phi)*sin(phi)^2)/r^4 - ...
              (R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/(2*r^4) + ...
              (15*R^3*cos(phi)^2*sin(phi)*(c31*cos(l) + s31*sin(l)))/r^4 - ...
              (6*R^2*cos(phi)*sin(phi)*(s22*sin(2*l) + c22*cos(2*l)))/r^3 + ...
              (3*R^2*c20*cos(phi)*sin(phi))/r^3);
 

    % FINAL DYNAMICS
    dr = v; 
    dl = csi; 
    dphi = n; 
    dv   = -mue/r^2 + r*n^2 + r*(csi + om)^2 *(cos(phi))^2 + ...
            aPg_r - ...
            l_v; 

    dcsi =  2*n*(csi + om)*tan(phi) - 2*v/r *(csi + om) + ...
            aPg_l/(r * cos(phi))^2 - ...
            l_csi/(r * cos(phi))^2; 

    dn   = -2*v/r *n - (csi + om)^2 *sin(phi)*cos(phi) + ...
            aPg_phi/r^2 - ...
            l_n/r^2; 

    dl_r = - l_n*((2*l_n)/r^3 + (2*n*v)/r^2 - (2*mue*((15*R^3*cos(phi)^3*(s32*sin(2*l) + c32*cos(2*l)))/r^4 + (3*R^2*cos(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 - (3*R^2*sin(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 + (R^3*c30*cos(phi)*(5*sin(phi)^2 - 3))/(2*r^4) - (30*R^3*cos(phi)*sin(phi)^2*(s32*sin(2*l) + c32*cos(2*l)))/r^4 - (45*R^3*cos(phi)^2*sin(phi)*(s33*sin(3*l) + c33*cos(3*l)))/r^4 + (5*R^3*c30*cos(phi)*sin(phi)^2)/r^4 - (R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/(2*r^4) + (15*R^3*cos(phi)^2*sin(phi)*(c31*cos(l) + s31*sin(l)))/r^4 - (6*R^2*cos(phi)*sin(phi)*(s22*sin(2*l) + c22*cos(2*l)))/r^3 + (3*R^2*c20*cos(phi)*sin(phi))/r^3))/r^3 - (mue*((60*R^3*cos(phi)^3*(s32*sin(2*l) + c32*cos(2*l)))/r^5 + (9*R^2*cos(phi)^2*(c21*cos(l) + s21*sin(l)))/r^4 - (9*R^2*sin(phi)^2*(c21*cos(l) + s21*sin(l)))/r^4 + (2*R^3*c30*cos(phi)*(5*sin(phi)^2 - 3))/r^5 - (120*R^3*cos(phi)*sin(phi)^2*(s32*sin(2*l) + c32*cos(2*l)))/r^5 - (180*R^3*cos(phi)^2*sin(phi)*(s33*sin(3*l) + c33*cos(3*l)))/r^5 + (20*R^3*c30*cos(phi)*sin(phi)^2)/r^5 - (2*R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/r^5 + (60*R^3*cos(phi)^2*sin(phi)*(c31*cos(l) + s31*sin(l)))/r^5 - (18*R^2*cos(phi)*sin(phi)*(s22*sin(2*l) + c22*cos(2*l)))/r^4 + (9*R^2*c20*cos(phi)*sin(phi))/r^4))/r^2) - l_v*(cos(phi)^2*(csi + om)^2 + mue*((6*R^2*c20*(3*sin(phi)^2 - 1))/r^5 + (36*R^2*cos(phi)^2*(s22*sin(2*l) + c22*cos(2*l)))/r^5 + (300*R^3*cos(phi)^3*(s33*sin(3*l) + c33*cos(3*l)))/r^6 + (10*R^3*c30*sin(phi)*(5*sin(phi)^2 - 3))/r^6 + (36*R^2*cos(phi)*sin(phi)*(c21*cos(l) + s21*sin(l)))/r^5 + (300*R^3*cos(phi)^2*sin(phi)*(s32*sin(2*l) + c32*cos(2*l)))/r^6 + (10*R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/r^6) + (2*mue)/r^3 + n^2) - l_csi*((2*l_csi)/(r^3*cos(phi)^2) + (2*v*(csi + om))/r^2 + (2*mue*((3*R^2*cos(phi)^2*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^3 + (15*R^3*cos(phi)^3*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^4 + (3*R^2*cos(phi)*sin(phi)*(c21*sin(l) - s21*cos(l)))/r^3 + (15*R^3*cos(phi)^2*sin(phi)*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + (R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/(2*r^4)))/(r^3*cos(phi)^2) + (mue*((9*R^2*cos(phi)^2*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^4 + (60*R^3*cos(phi)^3*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^5 + (9*R^2*cos(phi)*sin(phi)*(c21*sin(l) - s21*cos(l)))/r^4 + (60*R^3*cos(phi)^2*sin(phi)*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^5 + (2*R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/r^5))/(r^2*cos(phi)^2));

    dl_l = (l_csi*mue*((3*R^2*cos(phi)^2*(4*s22*sin(2*l) + 4*c22*cos(2*l)))/r^3 + (15*R^3*cos(phi)^3*(9*s33*sin(3*l) + 9*c33*cos(3*l)))/r^4 + (3*R^2*cos(phi)*sin(phi)*(c21*cos(l) + s21*sin(l)))/r^3 + (15*R^3*cos(phi)^2*sin(phi)*(4*s32*sin(2*l) + 4*c32*cos(2*l)))/r^4 + (R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/(2*r^4)))/(r^2*cos(phi)^2) - (l_n*mue*((3*R^2*sin(phi)^2*(c21*sin(l) - s21*cos(l)))/r^3 - (3*R^2*cos(phi)^2*(c21*sin(l) - s21*cos(l)))/r^3 - (15*R^3*cos(phi)^3*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + (30*R^3*cos(phi)*sin(phi)^2*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + (45*R^3*cos(phi)^2*sin(phi)*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^4 + (R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/(2*r^4) - (15*R^3*cos(phi)^2*sin(phi)*(c31*sin(l) - s31*cos(l)))/r^4 + (6*R^2*cos(phi)*sin(phi)*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^3))/r^2 - l_v*mue*((9*R^2*cos(phi)^2*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^4 + (60*R^3*cos(phi)^3*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^5 + (9*R^2*cos(phi)*sin(phi)*(c21*sin(l) - s21*cos(l)))/r^4 + (60*R^3*cos(phi)^2*sin(phi)*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^5 + (2*R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/r^5); 

    dl_phi = l_n*(cos(phi)^2*(csi + om)^2 - sin(phi)^2*(csi + om)^2 + (mue*((6*R^2*cos(phi)^2*(s22*sin(2*l) + c22*cos(2*l)))/r^3 + (45*R^3*cos(phi)^3*(s33*sin(3*l) + c33*cos(3*l)))/r^4 - (3*R^2*c20*cos(phi)^2)/r^3 - (6*R^2*sin(phi)^2*(s22*sin(2*l) + c22*cos(2*l)))/r^3 - (30*R^3*sin(phi)^3*(s32*sin(2*l) + c32*cos(2*l)))/r^4 + (3*R^2*c20*sin(phi)^2)/r^3 + (5*R^3*c30*sin(phi)^3)/r^4 - (15*R^3*cos(phi)^3*(c31*cos(l) + s31*sin(l)))/r^4 + (R^3*c30*sin(phi)*(5*sin(phi)^2 - 3))/(2*r^4) + (12*R^2*cos(phi)*sin(phi)*(c21*cos(l) + s21*sin(l)))/r^3 + (105*R^3*cos(phi)^2*sin(phi)*(s32*sin(2*l) + c32*cos(2*l)))/r^4 - (90*R^3*cos(phi)*sin(phi)^2*(s33*sin(3*l) + c33*cos(3*l)))/r^4 + (R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/(2*r^4) - (15*R^3*c30*cos(phi)^2*sin(phi))/r^4 + (45*R^3*cos(phi)*sin(phi)^2*(c31*cos(l) + s31*sin(l)))/r^4))/r^2) - l_csi*(2*n*(csi + om)*(tan(phi)^2 + 1) - (2*l_csi*sin(phi))/(r^2*cos(phi)^3) + (mue*((3*R^2*sin(phi)^2*(c21*sin(l) - s21*cos(l)))/r^3 - (3*R^2*cos(phi)^2*(c21*sin(l) - s21*cos(l)))/r^3 - (15*R^3*cos(phi)^3*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + (30*R^3*cos(phi)*sin(phi)^2*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + (45*R^3*cos(phi)^2*sin(phi)*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^4 + (R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/(2*r^4) - (15*R^3*cos(phi)^2*sin(phi)*(c31*sin(l) - s31*cos(l)))/r^4 + (6*R^2*cos(phi)*sin(phi)*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^3))/(r^2*cos(phi)^2) - (2*mue*sin(phi)*((3*R^2*cos(phi)^2*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^3 + (15*R^3*cos(phi)^3*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^4 + (3*R^2*cos(phi)*sin(phi)*(c21*sin(l) - s21*cos(l)))/r^3 + (15*R^3*cos(phi)^2*sin(phi)*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + (R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/(2*r^4)))/(r^2*cos(phi)^3)) + l_v*(mue*((60*R^3*cos(phi)^3*(s32*sin(2*l) + c32*cos(2*l)))/r^5 + (9*R^2*cos(phi)^2*(c21*cos(l) + s21*sin(l)))/r^4 - (9*R^2*sin(phi)^2*(c21*cos(l) + s21*sin(l)))/r^4 + (2*R^3*c30*cos(phi)*(5*sin(phi)^2 - 3))/r^5 - (120*R^3*cos(phi)*sin(phi)^2*(s32*sin(2*l) + c32*cos(2*l)))/r^5 - (180*R^3*cos(phi)^2*sin(phi)*(s33*sin(3*l) + c33*cos(3*l)))/r^5 + (20*R^3*c30*cos(phi)*sin(phi)^2)/r^5 - (2*R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/r^5 + (60*R^3*cos(phi)^2*sin(phi)*(c31*cos(l) + s31*sin(l)))/r^5 - (18*R^2*cos(phi)*sin(phi)*(s22*sin(2*l) + c22*cos(2*l)))/r^4 + (9*R^2*c20*cos(phi)*sin(phi))/r^4) + 2*r*cos(phi)*sin(phi)*(csi + om)^2); 
    
    dl_v = (2*l_n*n)/r - l_r + (2*l_csi*(csi + om))/r;

    dl_csi = l_csi*((2*v)/r - 2*n*tan(phi)) - l_l - l_v*r*cos(phi)^2*(2*csi + 2*om) + l_n*cos(phi)*sin(phi)*(2*csi + 2*om);

    dl_n = (2*l_n*v)/r - l_phi - 2*l_v*n*r - 2*l_csi*tan(phi)*(csi + om);
 
    
    % output
    dsph = [dr; dl; dphi; dv; dcsi; dn; dl_r; dl_l; dl_phi; dl_v; dl_csi; dl_n]; 

end



function [cost_dyn] = costate 
    % Compute costate dynamics
    syms r l phi v csi n mue R om l_r l_l l_phi l_v l_csi l_n c20 c21 c22 s33 c30 c31 c32 c33 s21 s22 s31 s32 s33
    
    aPg_r = -mue*((3*R^2*c20*(3*sin(phi)^2 - 1))/(2*r^4) + ...
            (9*R^2*cos(phi)^2*(s22*sin(2*l) + c22*cos(2*l)))/r^4 + ...
            (60*R^3*cos(phi)^3*(s33*sin(3*l) + c33*cos(3*l)))/r^5 + ...
            (2*R^3*c30*sin(phi)*(5*sin(phi)^2 - 3))/r^5 + ...
            (9*R^2*cos(phi)*sin(phi)*(c21*cos(l) + s21*sin(l)))/r^4 + ...
            (60*R^3*cos(phi)^2*sin(phi)*(s32*sin(2*l) + c32*cos(2*l)))/r^5 + ...
            (2*R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/r^5);
    
    
    
    
    aPg_l =  -mue*((3*R^2*cos(phi)^2*(2*c22*sin(2*l) - 2*s22*cos(2*l)))/r^3 + ...
             (15*R^3*cos(phi)^3*(3*c33*sin(3*l) - 3*s33*cos(3*l)))/r^4 + ...
             (3*R^2*cos(phi)*sin(phi)*(c21*sin(l) - s21*cos(l)))/r^3 + ...
             (15*R^3*cos(phi)^2*sin(phi)*(2*c32*sin(2*l) - 2*s32*cos(2*l)))/r^4 + ...
             (R^3*cos(phi)*(15*sin(phi)^2 - 3)*(c31*sin(l) - s31*cos(l)))/(2*r^4));
    
    
    
    aPg_phi = mue*((15*R^3*cos(phi)^3*(s32*sin(2*l) + c32*cos(2*l)))/r^4 + ...
              (3*R^2*cos(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 - ...
              (3*R^2*sin(phi)^2*(c21*cos(l) + s21*sin(l)))/r^3 + ...
              (R^3*c30*cos(phi)*(5*sin(phi)^2 - 3))/(2*r^4) - ...
              (30*R^3*cos(phi)*sin(phi)^2*(s32*sin(2*l) + c32*cos(2*l)))/r^4 - ...
              (45*R^3*cos(phi)^2*sin(phi)*(s33*sin(3*l) + c33*cos(3*l)))/r^4 + ...
              (5*R^3*c30*cos(phi)*sin(phi)^2)/r^4 - ...
              (R^3*sin(phi)*(15*sin(phi)^2 - 3)*(c31*cos(l) + s31*sin(l)))/(2*r^4) + ...
              (15*R^3*cos(phi)^2*sin(phi)*(c31*cos(l) + s31*sin(l)))/r^4 - ...
              (6*R^2*cos(phi)*sin(phi)*(s22*sin(2*l) + c22*cos(2*l)))/r^3 + ...
              (3*R^2*c20*cos(phi)*sin(phi))/r^3);
    
    
    
    dr = v; 
    dl = csi; 
    dphi = n; 
    
    dv   = -mue/r^2 + r*n^2 + r*(csi + om)^2 *(cos(phi))^2 + ...
            aPg_r - ...
            l_v;
    
    dcsi =  2*n*(csi + om)*tan(phi) - 2*v/r *(csi + om) + ...
            aPg_l/(r * cos(phi))^2 - ...
            l_csi/(r * cos(phi))^2;
    
    dn   = -2*v/r *n - (csi + om)^2 *sin(phi)*cos(phi) + ...
            aPg_phi/r^2 -...
            l_n/r^2;
    
    
    DfDx = jacobian([dr; dl; dphi; dv; dcsi; dn], [r, l, phi, v, csi, n]); 
    
    cost_dyn = - DfDx.' * [l_r; l_l; l_phi; l_v; l_csi; l_n];
    
end



function f = opt_solv(l, state, t_CONTROL, dat, grav)

    xf = state.final; x0 = state.initial; 
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    initial = [x0; l]; 
    [~, xx] = ode78(@(t,sph) SPHmotionGG_CONTROL(t, sph, dat, grav), t_CONTROL, initial, options);
    f = [xx(end,1)-xf(1); xx(end,2)-xf(2); xx(end,3)-xf(3);  xx(end,4)-xf(4); xx(end,5)-xf(5); xx(end,6)-xf(6)];

end

