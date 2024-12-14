%% MAE 425 Class Project
clc, clear,
mu = 398600.4418;   % [ km/s^2 ]
er = 6378.137;      % [ km ]

%% Parameter LankSat G15 Chaning Perifocal to Equatiorial 
a_sat = 31378;       % semi-major axis [ km ] 
e_sat = 0;           % eccentri [ degreee ]
i_sat = 10;          % Inclination [ degree ]
w_sat = 0;           % Arguemnt of perigee [ degree ]
omega_sat = 30;       % RAAN [degree]
True_sat = 0;           % true analomy [degree]

Apn_sat =A_PN(omega_sat,w_sat,i_sat);

v_sat = [-sqrt(mu/a_sat)*sind(True_sat);
    sqrt(mu/a_sat)*(e_sat+cosd(True_sat)); 
    0 ];

r_sat = [a_sat*cosd(True_sat);
    a_sat*sind(True_sat);
    0];

r_sat_1 = Apn_sat*r_sat;
v_sat_1 = Apn_sat*v_sat;
%% Parameters for G15 Chaning Perifocal to Equatiorial 
a_g = 42164.0;         % semi-major axis [ km ] 
e_g = 0;                % eccentri [ degreee ]
i_g = 0;                % Inclination [ degree ]
w_g = 0;
omega_g=0;
lamda_g = 90;            % true analomy [degree]

Apn_G15 = A_PN(omega_g,w_g,i_g);

v_15 = [-sqrt(mu/a_g)*sind(lamda_g);
    sqrt(mu/a_g)*(e_g+cosd(lamda_g)); 
    0 ];
r_G15 = [a_g*cosd(lamda_g);
    a_g*sind(lamda_g);
    0];

r_G15_1 = Apn_G15*r_G15;
v_G15_1 = Apn_G15*v_15;

%% Inlination change and Orbit Change

% tranferin orbit from LankSat orbit to G15 orbit
v_lsat = sqrt(mu/a_sat);                        %intial velocity at LankSat orbit
at = (a_sat + a_g)/2;                           % Finding the apogee between Lansat and G15
v_tranfer_sat = sqrt(mu *(2/a_sat - 1/at));     % velocity transfer orbit lanksat
delta_v_sat = abs( v_lsat - v_tranfer_sat);      % Change in velocity;

v_tranfer_G = sqrt(mu *(2/a_g - 1/at));         % velocity transfer orbit at G15
v_g = sqrt(mu/a_g);                             % orbit velocity at G15
delta_v_G15 = abs( v_tranfer_G -v_g);           % Velocity tranfer orbit of G15
delta_v_G15_i = 2*v_g*sind(i_sat/2) ;             % Velocity of (G-15) and v=2*v59*sin(I/2) inlination change

V_total = delta_v_sat + delta_v_G15 + delta_v_G15_i;   % Total Velocity Change
t_tranfer = pi*at^(3/2)/sqrt(mu);                      % Total Time 

direction = v_sat_1/norm(v_sat_1);
direction_2 = -v_G15_1/norm(v_G15_1);
v_new_2 = v_G15_1 + delta_v_G15*direction_2;
v_new = v_sat_1 + ((delta_v_sat)*direction);


dt = 1;
t=0;
nt = ceil(t_tranfer/dt);
tranfer = [r_sat_1;v_new];
tranfer_1 = rk4(tranfer,nt,t,dt,mu);

tranfer_G15 = [r_G15_1;v_new_2];
tranfer_G15_1 = rk4(tranfer_G15,nt,t,dt,mu);


%% Orbit intergation
% t=0;
% dt =1.0;
nt_ch = ceil(2*pi/sqrt(mu/a_sat^3)/dt);
x_chase=[r_sat_1;v_sat_1];
x_chase = rk4(x_chase,nt_ch,t,dt,mu);

dt = 1;
t=0;


nt_traget = ceil(2*pi/sqrt(mu/a_g^3)/dt);
x_traget(:,1)=[r_G15_1;v_G15_1];
x_traget = rk4(x_traget,nt_traget,t,dt,mu);

figure(1)
hold on
plot3(x_chase(1,:), x_chase(2,:), x_chase(3,:), 'g-', 'LineWidth',1);
plot3(x_traget(1,:), x_traget(2,:), x_traget(3,:), 'g-', 'LineWidth',1)
plot3(tranfer_1(1,:), tranfer_1(2,:), tranfer_1(3,:),'r', 'LineWidth',1)
%plot3(tranfer_G15_1(1,:),tranfer_G15_1(2,:),tranfer_G15_1(3,:),'k','LineWidth',1)
view(3)
%% Visualization






% Plot
% figure;
% hold on; grid on; axis equal;
% plot(r_circular_sat(1, :), r_circular_sat(2, :), 'b', 'LineWidth', 1.5); % LankSat orbit
% plot(r_circular_g15(1, :), r_circular_g15(2, :), 'g', 'LineWidth', 1.5); % G15 orbit
% plot(x_transfer, y_transfer, 'r--', 'LineWidth', 1.5);                  % Transfer orbit

% nt_traget = ceil(2*pi/sqrt(mu/a_g^3)/dt);
% x_traget(:,1)=[r_G15_1;v_G15_1];
% x_traget = rk4(nt_traget,x_traget,t,dt,mu);
% 
% 
% nt_ch_tranfer = ceil(t_tranfer/dt);
% for k = 1: nt_traget - 1
% 
%     % State at time K
%     x_k = xch(:,k);
% 
%     % RK 4 Step
%     k1 = twoBEOM(t, x_k, mu);
%     k2 = twoBEOM(t+dt/2, x_k + k1*dt/2, mu);
%     k3 = twoBEOM(t+dt/2, x_k + k2*dt/2, mu);
%     k4 = twoBEOM(t+dt, x_k + k3*dt, mu);
% 
%      x1_chase(:,k+1) = x_k + dt/6 * (k1+ 2*k2 + 2*k3 + k4);
% end
% nt2 = ceil(T_flight_vec(2)/dt);nt3 = ceil(2*T_flight_vec(3)/dt);
% nt4 = ceil(T_flight_vec(4)/dt);
% nt5 = ceil(2*T_flight_vec(5)/dt);nt6 = ceil(T_flight_vec(6)/dt);
% 
% x1(:,1) = [r_chase_1(:,1); v_chase_1(:,1)];
% x2(:,1) = [r_chase_1(:,2); v_chase_1(:,2)];
% x3(:,1) = [r_chase_1(:,3); v_chase_1(:,3)];
% x4(:,1) = [r_chase_1(:,4); v_chase_1(:,4)];
% x5(:,1) = [r_chase_1(:,5); v_chase_1(:,5)];
% x6(:,1) = [r_chase_1(:,6); v_chase_1(:,6)];
% 
% 
% 
% x_chase_1 = rk4(nt1,x1,t,dt,mu);
% x_chase_2 = rk4(nt2,x2,t,dt,mu);
% x_chase_3 = rk4(nt3,x3,t,dt,mu);
% x_chase_4 = rk4(nt4,x4,t,dt,mu);
% x_chase_5 = rk4(nt5,x5,t,dt,mu);
% x_chase_6 = rk4(nt6,x6,t,dt,mu);
% 
% 
% figure(1)
% 
% [xS, yS, zS] = sphere(50);
% 
% earth_radius = 6378137.0/1000; % m
% xSE = earth_radius*xS;
% ySE = earth_radius*yS;
% zSE = earth_radius*zS;
%  
% surface(xSE, ySE, zSE);
% 
% axis equal
% hold on
% plot3(x_chase(1,:), x_chase(2,:), x_chase(3,:), 'r-', 'LineWidth',1)
% plot3(x_traget(1,:), x_traget(2,:), x_traget(3,:), 'g-', 'LineWidth',1)
% plot3(x_chase_1(1,:), x_chase_1(2,:), x_chase_1(3,:), 'k-', 'LineWidth',1)
% plot3(x_chase_2(1,:), x_chase_2(2,:), x_chase_2(3,:), 'r-', 'LineWidth',1)
% plot3(x_chase_3(1,:), x_chase_3(2,:), x_chase_3(3,:), 'b-', 'LineWidth',1)
% plot3(x_chase_4(1,:), x_chase_4(2,:), x_chase_4(3,:), 'm-', 'LineWidth',1)
% plot3(x_chase_5(1,:), x_chase_5(2,:), x_chase_5(3,:), 'y-', 'LineWidth',1)
% plot3(x_chase_6(1,:), x_chase_6(2,:), x_chase_6(3,:), 'k-', 'LineWidth',1)
% 
% view(3)
% title('Lansat and G_15','FontSize',24)
% grid on
% hold off
% xlabel ('Inertial x (km)','FontSize',18)
% ylabel ('Inertial y (km)','FontSize',18)
% zlabel ('Inertial z (km)','FontSize',18)
% 
% 
% legend('Earth','Lansat','G15')
% save('project_propagation_save','x_chase_1','x_chase_2','x_chase_3','x_chase_4','x_chase_5','x_chase_6','dt')





