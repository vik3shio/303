clc; clear;

m = 1400; %kg
a = 1.14; %m
b = 1.33; %m
Cf = 25000; %N/rad
Cr = 21000; %N/rad
Iz = 2420; %kgm^2
u = 75/3.6; %km/hr

% Define constants for dx2/d2t = Adx/dt + Bdel
A = [-(Cf+Cr)/(m*u), -(a*Cf-b*Cr)/(m*u)-u;
       -0.0113, -((a^2)*Cf+(b^2)*Cr)/(Iz*u)];

del = 0.1;
B = [Cf/m; (a*Cf)/Iz];
B = del.*B;

% Compute y(t) and psi(t)

dt = 0.001; 
t = 0:dt:5;

x = zeros(4,length(t));
% dv/dt = d^2y/dt^2 = A(1,1)v + A(1,2)w + del*B
% dw/dt = d^2ψ/dt^2 = A(2,1)v + A(2,2)w + del*B

% IC at t = 0 (given eq7)
x(1,1) = 0; %y     
x(2,1) = 0; %psi  
x(3,1) = 0; %v 
x(4,1) = 0; %w

F = zeros(4,1);

F_temp = zeros(2,1);
xy_plot = zeros(2,length(t));

for n = 1:length(t)-1
    
    F_temp = [u*cos(x(2,n)) - x(3,n)+a*x(4,n)*sin(x(2,n));
                x(3,n)+a*x(4,n)*cos(x(2,n)) + u*sin(x(2,n))];

    F= [x(3,n);
        x(4,n);
        A(1,1)*x(3,n) + A(1,2)*x(4,n) + B(1);
        A(2,1)*x(3,n) + A(2,2)*x(4,n) + B(2)];

    x(:,n+1) = x(:,n) + dt * F(:);
    
    xy_plot(:,n+1) = xy_plot(:,n) + dt*F_temp(:);

end

figure;
hold on;
plot(xy_plot(1,:), xy_plot(2,:), 'b', 'LineWidth', 2); 
grid on;
xlabel('x');
ylabel('y');
