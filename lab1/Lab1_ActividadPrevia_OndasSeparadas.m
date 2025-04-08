clc; clear; close all;

A = 1;             
fc = 1000;        
Ts = 1/100000;     
t = 0:Ts:5/fc;    
m_t = A * sin(2 * pi * fc * t);  

fs = 5000;        
Ts_pam = 1/fs;     
d = 0.2;           
tau = d * Ts_pam;  

%% Tren de pulsos y PAM natural
pulsos_natural = zeros(size(t));  
m_pam_natural = zeros(size(t));   

for i = 1:length(t)
    if mod(t(i), Ts_pam) < tau
        pulsos_natural(i) = 1;  
        m_pam_natural(i) = m_t(i);
    end
end

%% Gráficas
figure;

subplot(3,1,1);
plot(t, m_t, 'k', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal Sinusoidal Original');
grid on;

subplot(3,1,2);
plot(t, pulsos_natural, 'r', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Tren de Pulsos (Muestreo Natural)');
grid on;

subplot(3,1,3);
plot(t, m_pam_natural, 'b', 'LineWidth', 1.5);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Señal PAM con Muestreo Natural (Pulsos Cuadrados)');
grid on;
