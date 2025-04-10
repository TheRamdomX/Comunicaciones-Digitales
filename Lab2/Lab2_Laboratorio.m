clc;
clear;
close all;

%% Parámetros de la señal
bit_rate = 1;                         % Tasa de bits (bps)
bits = randi([0 1], 1, 104);          % Secuencia aleatoria de bits
amplitud = 1;                         % Amplitud NRZ-L
muestras_por_bit = 40;                % Oversampling
Fs = bit_rate * muestras_por_bit;     % Frecuencia de muestreo
Ts = 1 / Fs;
t_total = length(bits) / bit_rate;
t = 0:Ts:t_total - Ts;

%% Codificación NRZ-L
senal_NRZ = repelem(2*bits - 1, muestras_por_bit); 

%% Parámetros del filtro de coseno alzado
roll_off_factors = [0, 0.25, 0.75, 1];
colores = ['b', 'r', 'g', 'm'];

for i = 1:length(roll_off_factors)
    alpha = roll_off_factors(i);

    % Filtro de coseno alzado
    span = 6;  % Duración del filtro en símbolos
    filtro = rcosdesign(alpha, span, muestras_por_bit, 'normal');

    % Filtrado
    senal_filtrada = conv(senal_NRZ, filtro, 'same');

    % Canal con ruido blanco
    SNR = 30;  % en dB
    senal_ruidosa = awgn(senal_filtrada, SNR, 'measured');

    % Diagrama de ojo
    eyediagram(senal_ruidosa, 2 * muestras_por_bit);
    title(['Diagrama de Ojo (α = ' num2str(alpha) ')']);
end
