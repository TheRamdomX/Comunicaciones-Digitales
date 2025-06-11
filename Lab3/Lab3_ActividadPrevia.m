%% Laboratorio 3 - Actividades Previas (ASK/OOK) - VERSIÓN FINAL CORREGIDA
clear all; close all; clc;

%% Parámetros de la señal
fc = 10000;       % Frecuencia portadora [Hz]
Rb = 1000;        % Tasa de bits [bps]
Tb = 1/Rb;        % Duración de un bit [s]
fs = 10*fc;       % Frecuencia de muestreo [Hz]
num_bits = 10;    % Número de bits a transmitir
t_total = num_bits*Tb; % Tiempo total de simulación

%% 1. Generación de señal modulante (tren de pulsos)
bits = randi([0 1], 1, num_bits);  % Secuencia aleatoria de bits
samples_per_bit = round(fs*Tb);
t = 0:1/fs:t_total-1/fs; % Vector de tiempo exacto

% Generar señal modulante con el tamaño exacto
m_t = zeros(1, length(t));
for i = 1:num_bits
    start_idx = (i-1)*samples_per_bit + 1;
    end_idx = i*samples_per_bit;
    if end_idx > length(t)
        end_idx = length(t);
    end
    m_t(start_idx:end_idx) = bits(i);
end

%% 2. Cálculo del ancho de banda teórico ASK
BW_ASK = 2*Rb;
fprintf('Ancho de banda teórico ASK: %.2f Hz\n', BW_ASK);

%% 3. Envolvente compleja para ASK
A = 1;  % Amplitud
g_t = A * m_t;  % Envolvente compleja

%% 4. Señal modulada ASK (OOK)
s_t = real(g_t .* exp(1j*2*pi*fc*t));

%% 5. Transformada de Fourier
N = length(t);
G_f = fftshift(fft(g_t));
f = (-N/2:N/2-1)*(fs/N);  % Vector de frecuencias

%% 6. Gráficos
figure;

subplot(3,1,1);
plot(t, m_t);
title('Señal Modulante m(t)'); xlabel('Tiempo [s]'); ylabel('Amplitud');
xlim([0 5*Tb]); grid on;

subplot(3,1,2);
plot(t, g_t);
title('Envolvente Compleja g(t) para ASK'); 
xlabel('Tiempo [s]'); ylabel('Amplitud');
xlim([0 5*Tb]); grid on;

subplot(3,1,3);
plot(f, abs(G_f)/max(abs(G_f)));
title('Espectro de la Envolvente Compleja G(f)');
xlabel('Frecuencia [Hz]'); ylabel('Magnitud Normalizada');
xlim([-3*Rb 3*Rb]); grid on;

%% Verificación ancho de banda
mask = abs(G_f) > 0.1*max(abs(G_f));
BW_medido = sum(mask) * (fs/N);
fprintf('Ancho de banda medido: %.2f Hz\n', BW_medido);

if abs(BW_medido - BW_ASK) < 0.2*BW_ASK
    disp('El ancho de banda teórico se cumple para ASK/OOK');
else
    disp('El ancho de banda teórico NO se cumple para ASK/OOK');
end


%% Para ASK/OOK - Versión mejorada de verificación
% Usar múltiples criterios de ancho de banda
BW_3dB = sum(abs(G_f) > 0.707*max(abs(G_f))) * (fs/N);
BW_10dB = sum(10*log10(abs(G_f)/max(abs(G_f))) > -10) * (fs/N);
BW_20dB = sum(10*log10(abs(G_f)/max(abs(G_f))) > -20) * (fs/N);

fprintf('Ancho de banda a -3 dB: %.2f Hz\n', BW_3dB);
fprintf('Ancho de banda a -10 dB: %.2f Hz\n', BW_10dB);
fprintf('Ancho de banda a -20 dB: %.2f Hz\n', BW_20dB);