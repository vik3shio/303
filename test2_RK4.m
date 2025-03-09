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

% IC at t = 0 (given eq7)
y0 = 0;
psi0 = 0;
v0 = -13.0964 + 24.4684 - 11.3720;     
w0 = -0.2496 - 0.6962 + 0.9457; 

%vector for i, i+1
x = zeros(4,length(t));
x(:,1) = [y0;psi0;v0;w0];
%vector for i+0.5 (intermediate steps)
xtemp = zeros(size(x));
%vector array for slopes
f = {};

%slope at i
f{1} = [x(3); x(4); A(1,1)*x(3) + A(1,2)*x(4) + B(1);
        A(2,1)*x(3) + A(2,2)*x(4) + B(2)];
    
for n = 1:length(t)-1

    %find xi+.5 and slope at i+.5
    xtemp = x(:,n) + 0.5*dt*f{1};
    f{2} = [xtemp(3); xtemp(4); 
        A(1,1)*xtemp(3) + A(1,2)*xtemp(4) + B(1);
        A(2,1)*xtemp(3) + A(2,2)*xtemp(4) + B(2)];

    %new i+0.5 and slope
    xtemp = x(:,n) + 0.5*dt*f{2};
    f{3} = [xtemp(3); xtemp(4); 
        A(1,1)*xtemp(3) + A(1,2)*xtemp(4) + B(1);
        A(2,1)*xtemp(3) + A(2,2)*xtemp(4) + B(2)];
    
    %find xi+1 and slope
    xtemp = x(:,n) + dt*f{3};
    f{4} = [xtemp(3); xtemp(4); 
        A(1,1)*xtemp(3) + A(1,2)*xtemp(4) + B(1);
        A(2,1)*xtemp(3) + A(2,2)*xtemp(4) + B(2)];
    
    f{5} = (1/6) .* f{1} + (1/3) .* f{2} + (1/3) .* f{3} + (1/2) * f{4};
    xtemp = x(:,n) + dt*f{5};

    x(:,n+1) = xtemp;
    f{1} = f{5};

end

y_exact = -13.0964 * exp(-1.9745 * t) + 24.4684 * exp(-0.9839 * t) - 11.3720;
psi_exact = -0.2496 * exp(-1.9745 * t) - 0.6962 * exp(-0.9839 * t) + 0.9457;

figure;
hold on;
plot(t, x(3,:), 'b', 'LineWidth', 2); 
plot(t, x(4,:), 'r', 'LineWidth', 2);
% plot(t,y_exact,'g', 'LineWidth', 2)
% plot(t,psi_exact,'g', 'LineWidth', 2)
grid on;
xlabel('Time (s)');
ylabel('Lateral Velocity and Yaw Rate');
legend('$\dot{y}(t)$','$\dot{\psi}(t)$', ...
    'eq(7) soln', 'Interpreter','LaTex');