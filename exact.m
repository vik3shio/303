clc; clear;

% Define time range
t = 0:0.001:5; % From 0 to 5 with step size 0.1

% Compute y(t) and ψ(t)
y = -13.0964 * exp(-1.9745 * t) + 24.4684 * exp(-0.9839 * t) - 11.3720;
psi = -0.2496 * exp(-1.9745 * t) - 0.6962 * exp(-0.9839 * t) + 0.9457;

% Plot y(t)
figure;
plot(t, y, 'b', 'LineWidth', 2);
hold on;
plot(t,psi)
grid on;
xlabel('Time (s)');
ylabel('y(t) and \psi(t)');
legend('y(t)', '\psi(t)');