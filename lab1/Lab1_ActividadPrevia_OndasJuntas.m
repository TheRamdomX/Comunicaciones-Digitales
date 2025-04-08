clear; clc; close all;

%% señal original
A = 1;                  % Amplitud 
fc = 1000;              % Frecuencia 
T_cont = 1/100000;      % Período de muestreo
t_final = 0.01;         % Tiempo final
t = 0:T_cont:t_final;   % Vector tiempo

m = A*sin(2*pi*fc*t);

%% PAM Natural
fs_natural = 5000;                       % Frecuencia
Ts_natural = 1/fs_natural;               % Período
d_natural = 0.8;                         % Ciclo de trabajo
tau_natural = d_natural * Ts_natural;    % Ancho del pulso

%% PAM Instantáneo
fs_instant = 5000;                       % Frecuencia 
Ts_instant = 1/fs_instant;               % Período 
d_instant = 0.1;                         % Ciclo de trabajo  (tau/T)
T_signal = 1/fc;                         % Período de la señal sinusoidal
tau_instant = d_instant * T_signal;      % Ancho 

%% Generación de la modulación PAM Natural
pam_natural = zeros(size(t));
t_sample_natural = 0:Ts_natural:t(end);
for k = 1:length(t_sample_natural)
    [~, idx] = min(abs(t - t_sample_natural(k)));
    pulse_indices = find(t >= t_sample_natural(k) & t < t_sample_natural(k) + tau_natural);
    pam_natural(pulse_indices) = m(idx);
end

%% Generación de la modulación PAM Instantáneo
pam_instant = zeros(size(t));
t_sample_instant = 0:Ts_instant:t(end);
for k = 1:length(t_sample_instant)
    [~, idx] = min(abs(t - t_sample_instant(k)));
    pulse_indices = find(t >= t_sample_instant(k) & t < t_sample_instant(k) + tau_instant);
    pam_instant(pulse_indices) = m(idx);
end

%% Graficar
figure;
plot(t, m, 'k', 'LineWidth', 1.5); hold on;
stem(t_sample_natural, A*sin(2*pi*fc*t_sample_natural), 'r', 'filled','LineWidth',1.5);
stairs(t, pam_instant, 'b', 'LineWidth', 1.5);
legend('Señal Original', 'PAM Natural', 'PAM Instantáneo');
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Modulación por Amplitud de Pulso (PAM)');
grid on;