clc; clear; close all;

%% Definición de la señal original
A = 1;                               % Amplitud de la señal senoidal
fc = 1000;                           % Frecuencia de la señal (Hz)
Ts = 1/100000;                       % Período de muestreo muy pequeño para simular señal continua
t = 0:Ts:5/fc;                       % Vector de tiempo para 5 períodos de la señal
m_t = A * sin(2 * pi * fc * t);      % Señal original senoidal

%% Modulación PAM Instantánea
fs = 5000;                      % Frecuencia de muestreo para PAM (Hz)
Ts_pam = 1/fs;                  % Período de muestreo PAM
m_pam_inst = zeros(size(t));    % Inicialización del vector PAM

% Selección de muestras en los instantes PAM
for i = 1:length(t)
    if mod(t(i), Ts_pam) < Ts
        m_pam_inst(i) = m_t(i);  % Se toma el valor de la señal en el instante PAM
    end
end

%% Cuantificación PCM
N = 2;                          % Número de bits de cuantificación
L = 2^N;                        % Número de niveles de cuantificación
m_max = max(abs(m_pam_inst));   % Valor máximo de la señal PAM
delta = 2 * m_max / L;          % Paso de cuantificación

% Cuantificación uniforme redondeando al nivel más cercano
m_pcm = delta * floor(m_pam_inst / delta + 0.5);

% Error de cuantificación
error_quant = m_pam_inst - m_pcm;

%% Visualización de resultados
% Subgráfico 1: Comparación entre señal original y PAM
subplot(3,1,1);
plot(t, m_t, 'k', t, m_pam_inst, 'b', 'LineWidth', 1.5);
title('Señal Original vs PAM Instantáneo');
legend('Original', 'PAM Instantáneo');
grid on;

% Subgráfico 2: Cuantificación PCM
subplot(3,1,2);
stairs(t, m_pcm, 'r', 'LineWidth', 1.5); % Señal PCM escalonada
hold on;
plot(t, m_pam_inst, 'k--');             % Señal PAM como referencia
title(['Señal Cuantificada (N = ', num2str(N), ' bits)']);
legend('Señal Cuantificada', 'PAM Instantáneo');
grid on;

% Subgráfico 3: Histograma del error de cuantificación
subplot(3,1,3);
histogram(error_quant, 10); 
title('Error de cuantificacion');
xlabel('Error');
ylabel('Frecuencia');
grid on;
