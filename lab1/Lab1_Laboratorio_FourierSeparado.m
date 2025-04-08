clc; clear; close all;

%% Parámetros de la señal original
A = 1;           
fc = 1000;         
Ts = 1/100000;     
t = 0:Ts:5/fc;     
m_t = A * sin(2 * pi * fc * t);

%% Parámetros de PAM
fs = 5000;         
Ts_pam = 1/fs;     
d = 0.2;            
tau = d * Ts_pam;  

%% Generación de la modulación PAM Natural
pulsos_natural = zeros(size(t));  
for i = 1:length(t)
    if mod(t(i), Ts_pam) < tau
        pulsos_natural(i) = 1; 
    end
end
m_pam_natural = m_t .* pulsos_natural; 

%% Generación de la modulación PAM Instantáneo
m_pam_inst = zeros(size(t));  
for i = 1:length(t)
    if mod(t(i), Ts_pam) < Ts  
        m_pam_inst(i) = m_t(i);     end
end

%% Transformadas de Fourier
N = length(t);                         
f = (0:N-1)*(1/Ts)/N;                   
f_pos = f(1:N/2);                       

M_t = abs(fft(m_t)/N);                  
M_pam_natural = abs(fft(m_pam_natural)/N);
M_pam_inst = abs(fft(m_pam_inst)/N);   

%% Graficar transformadas de Fourier
figure;
subplot(3,1,1);
plot(f_pos, M_t(1:N/2), 'k', 'LineWidth', 1.5);
title('FFT de la Señal Original');
xlabel('Frecuencia (Hz)'); ylabel('|M(f)|');
grid on;

subplot(3,1,2);
plot(f_pos, M_pam_natural(1:N/2), 'r', 'LineWidth', 1.5);
title('FFT de PAM Natural');
xlabel('Frecuencia (Hz)'); ylabel('|M(f)|');
grid on;

subplot(3,1,3);
plot(f_pos, M_pam_inst(1:N/2), 'b', 'LineWidth', 1.5);
title('FFT de PAM Instantáneo');
xlabel('Frecuencia (Hz)'); ylabel('|M(f)|');
grid on;