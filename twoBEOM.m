function xdot = twoBEOM(t, x, mu)
% Unpack x
x1 = x(1:3);
x2 = x(4:6);

% x1 Norm
x1n = sqrt(x1'*x1);

xdot  = [x2
        -mu/(x1n^3)*x1];
end