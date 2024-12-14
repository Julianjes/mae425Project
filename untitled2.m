% Parameters
mu = 398600; % Gravitational parameter (km^3/s^2)
r1 = 7000; % Initial orbit radius (km)
r2 = 15000; % Final orbit radius (km)

% Initial velocity for transfer (at periapsis of transfer orbit)
v_periapsis = sqrt(2*mu/r1 - mu/((r1+r2)/2)); % Vis-viva equation

% Time span for simulation
tspan = [0 10000]; % Adjust as needed for full orbit (s)

% Initial conditions
x0 = r1; % Starting at periapsis (km)
y0 = 0; % Along x-axis (km)
vx0 = 0; % Velocity along y-axis (km/s)
vy0 = v_periapsis; % Velocity at periapsis (km/s)

% Pack initial conditions
initial_conditions = [x0, y0, vx0, vy0];

% Solve using ode45
[t, sol] = ode45(@orbital_motion, tspan, initial_conditions);

% Extract position
x = sol(:, 1);
y = sol(:, 2);

% Plot initial and final orbits
theta = linspace(0, 2*pi, 500);
orbit1_x = r1 * cos(theta);
orbit1_y = r1 * sin(theta);
orbit2_x = r2 * cos(theta);
orbit2_y = r2 * sin(theta);

figure;
plot(orbit1_x, orbit1_y, 'b-', 'LineWidth', 1.5); hold on; % Initial orbit
plot(orbit2_x, orbit2_y, 'g-', 'LineWidth', 1.5); % Final orbit
plot(x, y, 'r--', 'LineWidth', 1.5); % Transfer trajectory

% Formatting
axis equal;
xlabel('x (km)');
ylabel('y (km)');
title('Hohmann Transfer Using Numerical Integration');
legend('Initial Orbit', 'Final Orbit', 'Transfer Orbit');
grid on;
hold off;

% Orbital motion function (at the end of the file)
function dydt = orbital_motion(t, state)
    % Extract state variables
    x = state(1);
    y = state(2);
    vx = state(3);
    vy = state(4);
    
    % Gravitational parameter
    mu = 398600; % km^3/s^2
    
    % Compute distance from the central body
    r = sqrt(x^2 + y^2); 
    
    % Compute accelerations
    ax = -mu * x / r^3;
    ay = -mu * y / r^3;
    
    % Return derivatives
    dydt = [vx; vy; ax; ay];
end
 