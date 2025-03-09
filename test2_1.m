clc; clear;

m = 1400; %kg
a = 1.14; %m
b = 1.33; %m
Cf = 25000; %N/rad
Cr = 21000; %N/rad
Iz = 2420; %kgm^2
u = 75; %km/h

% Define constants for dx2/d2t = Adx/dt + Bdel
% A = [-(Cf+Cr)/(m*u), -(a*Cf+b*Cr)/(m*u)-75;
%    -(a*Cf+b*Cr)/(Iz*u), -((a^2)*Cf+(b^2)*Cr)/(Iz*u)];

A = [-1.5771, -20.8529;
        -0.0113, -1.3812];


%B = [Cf/m; a*Cf/Iz];
del = 1;

B = [17.8571; 11.7769];

% Compute y(t) and psi(t)

dt = 0.001; 
t = 0:dt:5;

y = zeros(size(t));
psi = zeros(size(t));
v = zeros(size(t)); % v = dy/dt
w = zeros(size(t)); % w = dpsi/dt
% dv/dt = d^2y/dt^2 = A(1,1)v + A(1,2)w + del*B
% dw/dt = d^2ψ/dt^2 = A(2,1)v + A(2,2)w + del*B

% IC at t = 0 (given eq7)
y(1) = 0;     
psi(1) = 0;   
v(1) = -13.0964 + 24.4684 - 11.3720;     
w(1) = -0.2496 - 0.6962 + 0.9457; 

F = zeros(4,1);

for n = 1:length(t)-1
    
    F= [v(n);
        w(n);
        A(1,1)*v(n) + A(1,2)*w(n) + del*B(1,1);
        A(2,1)*v(n) + A(2,2)*w(n) + del*B(2,1)];

    y(n+1) = y(n) + dt * F(1,1);
    psi(n+1) = psi(n) + dt * F(2,1);
    v(n+1) = v(n) + dt * F(3,1);
    w(n+1) = w(n) + dt * F(4,1);

end

figure;
plot(t, v, 'b', 'LineWidth', 2); 
hold on;
plot(t, w, 'r', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('y(t) and \psi(t) Derivatives');
legend('$\dot{y}(t)$','$\dot{\psi}(t)$', 'Interpreter','LaTex');