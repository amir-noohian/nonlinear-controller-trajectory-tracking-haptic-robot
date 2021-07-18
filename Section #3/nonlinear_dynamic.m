function dx = nonlinear_dynamic(t,x)

global u1 u2 u3

dx = x;

dx(1:3) = x(4:6);

theta1 = x(1); theta2 = x(2); theta3 = x(3);
thetad1 = x(4); thetad2 = x(5); thetad3 = x(6);

% theta vector
thetad = [thetad1; thetad2; thetad3];

% U vector (input torques)
U = [u1; u2; u3];

%% model: the second model has been used

% the value of parameters
ma = 17.5e-3;
mc = 10.4e-3;
mbe = 0.2214;
mdf = 0.1106;
L1 = 13.97e-2;
L2 = 13.97e-2;
L3 = 0.325e-2;
L4 = 0.368e-2;
L5 = 0.527e-2;
g = 9.81;

% the elements of M matrix
m11 = ((0.5 * L1^2 + 0.125 * L2^2) * ma + (0.125 * L1^2 + 0.5 * L3^2) * mc) + 0.125 * L1^2 * (4*ma + mc) * cos(2*theta2) - 0.125 * (L2^2 * ma + 4 * L3^2 * mc) * cos(2*theta2) + 0.125 * (L2 * ma + L3 * mc) * cos(theta2) * sin(theta3);
m12 = 0;
m13 = 0;
m21 = m12;
m22 = L1^2 * (ma + 0.25 * mc);
m23 = -0.5 * L1 * (L2 * ma + L3 * mc) * sin(theta2 - theta3);
m31 = m13;
m32 = m23;
m33 = 0.25 * L2^2 * ma + L3^2 * mc;
% M matrix
M = [m11 m12 m13; m21 m22 m23; m31 m32 m33];

% the elements of V matrix
v11_1 = 0.25 * cos(theta3) * (2 * L1 * (L2 * ma + L3 * mc) * cos(theta2) + (L2^2 * ma + 4 * L3^2 * mc) * sin(theta3)) * thetad3;
v11 = 0.25 * (-2 * sin(theta2) * (L1^2 * (4 * ma + mc) * cos(theta2) + 2 * L1 * (L2 * ma + L3 * mc) * sin(theta3)) * thetad2) + v11_1;
v12 = -0.25 * (L1^2 * (4 * ma + mc) * sin(2 * theta2) + 2 * L1 * (L2 * ma + L3 * mc) * sin(theta2) * sin(theta3)) * thetad1;
v13 = -0.125 * (-4 * L1 * (L2 * ma + L3 * mc) * cos(theta2) * cos(theta3) - (L2^2 * ma + 4 * L3^2 * mc) * sin(2 * theta3)) * thetad1;
v21 = -v12;
v22 = 0;
v23 = 0.5 * L1 * (L2 * ma + L3 * mc) * cos(theta2 - theta3) * thetad3;
v31 = -v13;
v32 = 0.5 * L1 * (L2 * ma + L3 * mc) * cos(theta2 - theta3) * thetad2;
v33 = 0;
% V matrix
V = [v11 v12 v13; v21 v22 v23; v31 v32 v33];

% the elements of N vector
n1 = 0;
n2 = g * (L1 * (ma + 0.5 * mc) + L5 * mbe) * cos(theta2);
n3 = g * (0.5 * L2 * ma + L3 * mc - L4 * mdf) * sin(theta3);
% N vector
N = [n1; n2; n3];

%% 
B = V * thetad + N - U;

dx(4:6) = - M \ B;

