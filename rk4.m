function [xch] = rk4(nt,xch,t,dt,mu)

for k = 1: nt - 1

    % State at time K
    x_k = xch(:,k);

    % RK 4 Step
    k1 = twoBEOM(t, x_k, mu);
    k2 = twoBEOM(t+dt/2, x_k + k1*dt/2, mu);
    k3 = twoBEOM(t+dt/2, x_k + k2*dt/2, mu);
    k4 = twoBEOM(t+dt, x_k + k3*dt, mu);

     xch(:,k+1) = x_k + dt/6 * (k1+ 2*k2 + 2*k3 + k4);

    % Store new Value of x
   
end