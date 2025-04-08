clear; clc; close all;

%% Parámetros de la señal original
A = 1;
fc = 1000;
T_cont = 1/100000;
t_final = 0.01;
t = 0:T_cont:t_final;
m = A*sin(2*pi*fc*t);

%% Parámetros de PAM Natural
fs_natural = 5000;
Ts_natural = 1/fs_natural;
d_natural = 0.8;
tau_natural = d_natural * Ts_natural;

%% Parámetros de PAM Instantáneo
fs_instant = 5000;
Ts_instant = 1/fs_instant;
d_instant = 0.1;
T_signal = 1/fc;
tau_instant = d_instant * T_signal;

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

%% Transformadas de Fourier
N = length(t);
frecuencia = linspace(0, 1/(2*T_cont), N/2);

M_f = abs(fft(m)/N);
M_f = M_f(1:N/2);

PAM_N_f = abs(fft(pam_natural)/N);
PAM_N_f = PAM_N_f(1:N/2);

PAM_I_f = abs(fft(pam_instant)/N);
PAM_I_f = PAM_I_f(1:N/2);

%% Graficar transformadas de Fourier
figure;
plot(frecuencia, M_f, 'k', 'LineWidth', 1.5); hold on;
plot(frecuencia, PAM_N_f, 'r', 'LineWidth', 1.5);
plot(frecuencia, PAM_I_f, 'b', 'LineWidth', 1.5);
legend('Señal Original', 'PAM Natural', 'PAM Instantáneo');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
title('Transformada de Fourier de las Señales');
grid on;
