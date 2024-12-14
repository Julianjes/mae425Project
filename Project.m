%% MAE 425 Class Project
clc, clear,
mu = 398600.4418;   % [ km/s^2 ]
er = 6378.137;      % [ km ]

%% Parameter LankSat
LankSat.a = 31378;       % semi-major axis [ km ] 
LankSat.ecc = 0;           % eccentri [ degreee ]
LankSat.i = 10;          % Inclination [ degree ]
LankSat.w = 0;           % Arguemnt of perigee [ degree ]
LankSat.RAAN = 30;       % RAAN [degree]
LankSat.trueAn = 0;           % true analomy [degree]

Apn_sat =A_PN(LankSat);

v_sat= [-sqrt(mu/LankSat.a)*sind(LankSat.trueAn);
    sqrt(mu/LankSat.a)*(LankSat.ecc+cosd(LankSat.trueAn)); 
    0 ]; % Perifocal 

r_sat = [LankSat.a*cosd(LankSat.trueAn);
    LankSat.a*sind(LankSat.trueAn);
    0]; % Perifocal

r_satN = Apn_sat*r_sat; % Inertial
v_satN = Apn_sat*v_sat; % Inertial
%% Parameters for G15
G15.a = 42164.0;         % semi-major axis [ km ] 
G15.ecc = 0;                % eccentri [ degreee ]
G15.i = 0;                % Inclination [ degree ]
G15.w = 0;
G15.RAAN=0;
G15.truelon = 90;            % true analomy [degree]

Apn_G15 = A_PN(G15);

v_G15 = [-sqrt(mu/G15.a)*sind(G15.truelon);
    sqrt(mu/G15.a)*(G15.ecc+cosd(G15.truelon)); 
    0 ]; % Perifocal
r_G15 = [G15.a*cosd(G15.truelon);
    G15.a*sind(G15.truelon);
    0]; % Perifocal

r_G15N = Apn_G15*r_G15; % Inertial
v_G15N = Apn_G15*v_G15; % Inertial
%% Plane Change LankSat to G15, 10 degree inclination change

LankSat.V0 = sqrt(mu/LankSat.a); % LankSat Initial Velocity
deltai = 10; % [degrees]
dV_plane = 2 * v_satN * sind(deltai/2);

%% Hohmann Transfer Initial: a = 31378, desired: a = 42164
aT = (LankSat.a + G15.a)/2; % Transfer Apogee

vLankSat = sqrt(mu/LankSat.a); 
v_transferLank = sqrt(mu *(2/LankSat.a - 1/aT));
dV_transfer_Lanksat = abs(v_transferLank - vLankSat);

vG15 = sqrt(mu/G15.a); % Target Orbit
v_transferG15 = sqrt(mu *(2/G15.a- 1/aT));
dV_transfer_G15 = abs(vG15 - v_transferG15);

t_transfer = pi*aT^(3/2)/sqrt(mu);  % Transfer Period 

dV_plane_phase = dV_transfer_Lanksat + dV_transfer_G15 + dV_plane; % Total delta V for plane and transfer manuevers

%% Phase Manuever 

% 0.17 degrees behind G15

wtgt = sqrt(mu/(G15.a^3)); % angular rate of target
Ptgt = 2*pi/wtgt; % 
theta = 0.17;
k = 1;

dT = theta*pi/180/wtgt * k;
Pchase = Ptgt + dT; 
TOF = k * Pchase; % Time of Flight
a_chase = ((mu*Pchase^2)/(4*pi^2))^(1/3);       % Semi major axis 
v_phasing = sqrt(2/G15.a - 1/a_chase);            % Velocity of chase
dV_phase = abs(vG15 - v_phasing);           % Total velocity change in phase 

% Total delta v to perform all manuevers 
dV_total = dV_plane_phase + dV_phase;


%% Orbit intergation
t=0;
dt =1.0;
nt_ch = ceil(2*pi/sqrt(mu/LankSat.a^3)/dt);
x_chase = zeros(6,nt_ch);
x_chase(:,1)=[r_satN;v_satN];
x_chase = rk4(nt_ch,x_chase,t,dt,mu);

nt_target = ceil(2*pi/sqrt(mu/G15.a^3)/dt);
x_target(:,1)=[r_G15N;v_G15N];
x_target = rk4(nt_target,x_target,t,dt,mu);

nt_ch_tranfer = ceil(t_transfer/dt);

for k = 1: nt_target - 1

    % State at time K
    x_k = xch(:,k);

    % RK 4 Step
    k1 = twoBEOM(t, x_k, mu);
    k2 = twoBEOM(t+dt/2, x_k + k1*dt/2, mu);
    k3 = twoBEOM(t+dt/2, x_k + k2*dt/2, mu);
    k4 = twoBEOM(t+dt, x_k + k3*dt, mu);

    x1_chase(:,k+1) = x_k + dt/6 * (k1+ 2*k2 + 2*k3 + k4);
end
nt2 = ceil(T_flight_vec(2)/dt);nt3 = ceil(2*T_flight_vec(3)/dt);
nt4 = ceil(T_flight_vec(4)/dt);
nt5 = ceil(2*T_flight_vec(5)/dt);nt6 = ceil(T_flight_vec(6)/dt);

x1(:,1) = [r_chase_1(:,1); v_chase_1(:,1)];
x2(:,1) = [r_chase_1(:,2); v_chase_1(:,2)];
x3(:,1) = [r_chase_1(:,3); v_chase_1(:,3)];
x4(:,1) = [r_chase_1(:,4); v_chase_1(:,4)];
x5(:,1) = [r_chase_1(:,5); v_chase_1(:,5)];
x6(:,1) = [r_chase_1(:,6); v_chase_1(:,6)];



x_chase_1 = rk4(nt1,x1,t,dt,mu);
x_chase_2 = rk4(nt2,x2,t,dt,mu);
x_chase_3 = rk4(nt3,x3,t,dt,mu);
x_chase_4 = rk4(nt4,x4,t,dt,mu);
x_chase_5 = rk4(nt5,x5,t,dt,mu);
x_chase_6 = rk4(nt6,x6,t,dt,mu);


figure(1)

[xS, yS, zS] = sphere(50);

earth_radius = 6378137.0/1000; % m
xSE = earth_radius*xS;
ySE = earth_radius*yS;
zSE = earth_radius*zS;
 
surface(xSE, ySE, zSE);

axis equal
hold on
plot3(x_chase(1,:), x_chase(2,:), x_chase(3,:), 'r-', 'LineWidth',1)
plot3(x_target(1,:), x_target(2,:), x_target(3,:), 'g-', 'LineWidth',1)
plot3(x_chase_1(1,:), x_chase_1(2,:), x_chase_1(3,:), 'k-', 'LineWidth',1)
plot3(x_chase_2(1,:), x_chase_2(2,:), x_chase_2(3,:), 'r-', 'LineWidth',1)
plot3(x_chase_3(1,:), x_chase_3(2,:), x_chase_3(3,:), 'b-', 'LineWidth',1)
plot3(x_chase_4(1,:), x_chase_4(2,:), x_chase_4(3,:), 'm-', 'LineWidth',1)
plot3(x_chase_5(1,:), x_chase_5(2,:), x_chase_5(3,:), 'y-', 'LineWidth',1)
plot3(x_chase_6(1,:), x_chase_6(2,:), x_chase_6(3,:), 'k-', 'LineWidth',1)

view(3)
title('Lansat and G_15','FontSize',24)
grid on
hold off
xlabel ('Inertial x (km)','FontSize',18)
ylabel ('Inertial y (km)','FontSize',18)
zlabel ('Inertial z (km)','FontSize',18)


legend('Earth','Lansat','G15')
save('project_propagation_save','x_chase_1','x_chase_2','x_chase_3','x_chase_4','x_chase_5','x_chase_6','dt')





