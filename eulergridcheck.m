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

F = zeros(4,1);

dt = [1,0.5,0.1,0.01,0.001];
colors = ['r','b','g','c','k'];

for i =1:length(dt)

    t = 0:dt(i):5;
    
    x = zeros(4,length(t));
    % dv/dt = d^2y/dt^2 = A(1,1)v + A(1,2)w + del*B
    % dw/dt = d^2ψ/dt^2 = A(2,1)v + A(2,2)w + del*B
    
    % IC at t = 0 (given eq7)
    x(1,1) = 0; %y     
    x(2,1) = 0; %psi  
    x(3,1) = -13.0964 + 24.4684 - 11.3720; %v 
    x(4,1) = -0.2496 - 0.6962 + 0.9457; %w

    for n = 1:length(t)-1
        
        F= [x(3,n);
            x(4,n);
            A(1,1)*x(3,n) + A(1,2)*x(4,n) + B(1);
            A(2,1)*x(3,n) + A(2,2)*x(4,n) + B(2)];
    
        x(:,n+1) = x(:,n) + dt(i) * F(:,1);

    end

    subplot(2,1,1);
    plot(t, x(3,:), 'color',colors(i), 'LineWidth', 1); 
    xlabel('Time (s)');
    ylabel('Lateral Velocity');
    title('Lateral Velocity Grid Independence Check');
    legend('dt = 1', 'dt = 0.5', 'dt = 0.1', 'dt = 0.01','dt = 0.001');
    grid on;
    hold on;
    
    subplot(2,1,2);
    plot(t, x(4,:), 'color',colors(i), 'LineWidth', 1); 
    xlabel('Time (s)');
    ylabel('Yaw Rate');
    title('Yaw Rate Grid Independence Check');
    grid on;
    hold on;

end