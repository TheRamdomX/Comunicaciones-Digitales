clear all; close all; clc;

%% Parámetros de la señal
fc = 10000;       % Frecuencia portadora central [Hz]
delta_f = 2000;   % Desviación de frecuencia [Hz]
Rb = 1000;        % Tasa de bits [bps]
Tb = 1/Rb;        % Duración de un bit [s]
fs = 10*(fc+delta_f); % Frecuencia de muestreo [Hz]
num_bits = 10;    % Número de bits a transmitir
t_total = num_bits*Tb; % Tiempo total de simulación

%% 1. Generación de señal modulante
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

%% 2. Cálculo del ancho de banda teórico FSK
BW_FSK = 2*delta_f + 2*Rb;
fprintf('Ancho de banda teórico FSK: %.2f Hz\n', BW_FSK);

%% 3. Envolvente compleja para FSK
A = 1;  % Amplitud
int_m_t = cumsum(m_t)/fs; % Integral de m(t)
g_t = A * exp(1j*2*pi*delta_f*int_m_t);

%% 4. Señal modulada FSK
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
plot(t, real(g_t), 'b', t, imag(g_t), 'r');
title('Envolvente Compleja g(t) para FSK');
xlabel('Tiempo [s]'); ylabel('Amplitud');
xlim([0 5*Tb]); legend('Real','Imag'); grid on;

subplot(3,1,3);
plot(f, abs(G_f)/max(abs(G_f)));
title('Espectro de la Envolvente Compleja G(f)');
xlabel('Frecuencia [Hz]'); ylabel('Magnitud Normalizada');
xlim([-2*(delta_f+Rb) 2*(delta_f+Rb)]); grid on;

%% Análisis de espectro de la señal modulada
S_f = fftshift(fft(s_t));

figure;
plot(f, 10*log10(abs(S_f)/max(abs(S_f))));
title('Espectro de la señal FSK modulada (dB)');
xlabel('Frecuencia [Hz]'); ylabel('Magnitud [dB]');
xlim([fc-2*BW_FSK fc+2*BW_FSK]); grid on;

% Medición del ancho de banda
mask = 10*log10(abs(S_f)/max(abs(S_f))) > -20;
BW_medido = sum(mask) * (fs/N);
fprintf('Ancho de banda medido (-20 dB): %.2f Hz\n', BW_medido);

if abs(BW_medido - BW_FSK) < 0.2*BW_FSK
    disp('El ancho de banda teórico se cumple para FSK');
else
    disp('El ancho de banda teórico NO se cumple para FSK');
end