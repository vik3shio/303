clc; clear;

m = 1400; %kg
a = 1.14; %m
b = 1.33; %m
Cf = 25000; %N/rad
Cr = 21000; %N/rad
Iz = 2420; %kgm^2

del = 0.1;
B = [Cf/m; (a*Cf)/Iz];
B = del.*B;

% Compute y(t) and psi(t)
dt = 0.01; 
t = 0:dt:5;

x = zeros(4,length(t));

% IC at t = 0 (given eq7)
x(1,1) = 0; %y     
x(2,1) = 0; %psi  
x(3,1) = -13.0964 + 24.4684 - 11.3720; %v 
x(4,1) = -0.2496 - 0.6962 + 0.9457; %w
F = zeros(4,1);

u_var = [20,50,75,100];
colors = ['r','b','g','c','k','m',"#EDB120"];

figure
hold on;

for i = 1:length(u_var)
    u = u_var(i)/3.6; %m/s
    
    % Define constants for dx2/d2t = Adx/dt + Bdel
    A = [-(Cf+Cr)/(m*u), -(a*Cf-b*Cr)/(m*u)-u;
           -0.0113, -((a^2)*Cf+(b^2)*Cr)/(Iz*u)];
    
    for n = 1:length(t)-1
        
        F= [x(3,n);
            x(4,n);
            A(1,1)*x(3,n) + A(1,2)*x(4,n) + B(1);
            A(2,1)*x(3,n) + A(2,2)*x(4,n) + B(2)];
    
        x(:,n+1) = x(:,n) + dt * F(:,1);
        
        
    end
    
    subplot(2,1,1); 
    plot(t, x(3,:),'color', colors(i), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Lateral Accel. (m/s^2)');
    title(['Variable Tangential Velocity']);
    legend('u = 20 km/h', 'u = 50 km/h', 'u = 75 km/h', ...
        'u = 100 km/h','u = 200','u = 300');
    hold on;
    grid on;
    
    subplot(2,1,2);   
    plot(t, x(2,:),'color', colors(i),'LineWidth', 2); 
    xlabel('Time (s)');
    ylabel('Yaw Rate (rad/s)');
    hold on;
    grid on;

    %legend('y_a(t)','\psi_v(t)');
    x = zeros;
    x(1,1) = 0; %y     
    x(2,1) = 0; %psi  
    x(3,1) = -13.0964 + 24.4684 - 11.3720; %v 
    x(4,1) = -0.2496 - 0.6962 + 0.9457; %w

end

