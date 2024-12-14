%% Ground Track

mu = 398600.4418; % gravitational constant [km^3/s^2]
R_earth = 6378.137; % [km]

load('project_propagation_save.mat')
year= 2024;
month = 10;
day= 1;
hour = 0;
min=0;
sec=0;

JD = 367*(year) - floor(7*(year + floor((month+9)/12))/4)...
    + floor(275*month/9) + day + 1721013.5 + (((sec/60) + min)/60 + hour)/24; 

%% Generating Groud Trace

% Initializing, Lat, longitude, altitude, r_ECEF
nt1 = length(x_chase_1); nt2 = length(x_chase_2); nt3 = length(x_chase_3); 
nt4 = length(x_chase_4); nt5 = length(x_chase_5); nt6 = length(x_chase_6);

t_1 = 0:dt:nt1; t_2 = 0:dt:nt2; t_3 = 0:dt:nt3; t_4 = 0:dt:nt4;
t_5 = 0:dt:nt5; t_6 = 0:dt:nt6;

lat1 = zeros(1,nt1); lat2 = zeros(1,nt2); lat3 = zeros(1,nt3); lat4 = zeros(1,nt4);
lat5 = zeros(1,nt5); lat6 = zeros(1,nt6);

lon1 = zeros(1,nt1); lon2 = zeros(1,nt2); lon3 = zeros(1,nt3); lon4 = zeros(1,nt4);
lon5 = zeros(1,nt5); lon6 = zeros(1,nt6);

r_ECEF_1 = zeros(3,nt1); r_ECEF_2 = zeros(3,nt2); r_ECEF_3 = zeros(3,nt3); r_ECEF_4 = zeros(3,nt4);
r_ECEF_5 = zeros(3,nt5); r_ECEF_6 = zeros(3,nt6);

alt_1 = zeros(1,nt1); alt_2 = zeros(1,nt2); alt_3 = zeros(1,nt3); alt_4 = zeros(1,nt4);
alt_5 = zeros(1,nt5); alt_6 = zeros(1,nt6);

%% Looping
for i = 1: nt1
    ti = t_1(i);
    r_ECI_i = x_chase_1(1:3,i);
    JD_1 = JD + ti/86400;
    r_ECEF_1(1:3,i) = ECI_ECEF(JD_1,r_ECI_i);
    r_i = norm(r_ECEF_1(1:3,i));

    lat1(i) = asind(r_ECEF_1(3,i)/r_i);
    lon1(i) = atan2d(r_ECEF_1(2,i),r_ECEF_1(1,i));
    alt_1(i) = r_i - R_earth;

end
JD_11 = JD_1;

for i = 1: nt2
    ti = t_2(i);
    r_ECI_i = x_chase_2(1:3,i);
    JD_2 = JD_11 + ti/86400;
    r_ECEF_2(1:3,i) = ECI_ECEF(JD_2,r_ECI_i);
    r_i = norm(r_ECEF_2(1:3,i));

    lat2(i) = asind(r_ECEF_2(3,i)/r_i);
    lon2(i) = atan2d(r_ECEF_2(2,i),r_ECEF_1(1,i));
    alt_2(i) = r_i - R_earth;

end
JD_22 = JD_2;

for i = 1: nt3
    ti = t_3(i);
    r_ECI_i = x_chase_3(1:3,i);
    JD_3 = JD_22 + ti/86400;
    r_ECEF_3(1:3,i) = ECI_ECEF(JD_3,r_ECI_i);
    r_i = norm(r_ECEF_3(1:3,i));

    lat3(i) = asind(r_ECEF_3(3,i)/r_i);
    lon3(i) = atan2d(r_ECEF_3(2,i),r_ECEF_3(1,i));
    alt_3(i) = r_i - R_earth;

end

JD_33 = JD_3;

for i = 1: nt4
    ti = t_4(i);
    r_ECI_i = x_chase_4(1:3,i);
    JD_4 = JD_33 + ti/86400;
    r_ECEF_4(1:3,i) = ECI_ECEF(JD_4,r_ECI_i);
    r_i = norm(r_ECEF_4(1:3,i));

    lat4(i) = asind(r_ECEF_4(3,i)/r_i);
    lon4(i) = atan2d(r_ECEF_4(2,i),r_ECEF_4(1,i));
    alt_4(i) = r_i - R_earth;

end
JD_44 = JD_4;

for i = 1: nt5
    ti = t_5(i);
    r_ECI_i = x_chase_5(1:3,i);
    JD_5 = JD_44 + ti/86400;
    r_ECEF_5(1:3,i) = ECI_ECEF(JD_5,r_ECI_i);
    r_i = norm(r_ECEF_5(1:3,i));

    lat5(i) = asind(r_ECEF_5(3,i)/r_i);
    lon5(i) = atan2d(r_ECEF_5(2,i),r_ECEF_5(1,i));
    alt_5(i) = r_i - R_earth;

end

JD_55 = JD_5;

for i = 1: nt6
    ti = t_6(i);
    r_ECI_i = x_chase_6(1:3,i);
    JD_6 = JD_55 + ti/86400;
    r_ECEF_6(1:3,i) = ECI_ECEF(JD_6,r_ECI_i);
    r_i = norm(r_ECEF_6(1:3,i));

    lat6(i) = asind(r_ECEF_6(3,i)/r_i);
    lon6(i) = atan2d(r_ECEF_6(2,i),r_ECEF_6(1,i));
    alt_6(i) = r_i - R_earth;

end



%% Plotthing Ground Trace
figure(1)
hold on
grid on
plot(lon1,lat1,'.')
plot(lon2,lat2,'.')
plot(lon3,lat3,'.')
plot(lon4,lat4,'.')
plot(lon5,lat5,'.')
plot(lon6,lat6,'.')