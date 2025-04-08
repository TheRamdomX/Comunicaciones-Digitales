clc; clear; close all;


A = 1;              
fc = 1000;         
Ts = 1/100000;     
t = 0:Ts:5/fc;     
m_t = A * sin(2 * pi * fc * t); 

fs = 5000;          
Ts_pam = 1/fs;      
m_pam_inst = zeros(size(t)); 
for i = 1:length(t)
    if mod(t(i), Ts_pam) < Ts
        m_pam_inst(i) = m_t(i); 
    end
end

N = 2;              
L = 2^N;            
m_max = max(abs(m_pam_inst)); 
delta = 2 * m_max / L;

m_pcm = delta * floor(m_pam_inst / delta + 0.5);


error_quant = m_pam_inst - m_pcm;

figure;
subplot(3,1,1);
plot(t, m_t, 'k', t, m_pam_inst, 'b', 'LineWidth', 1.5);
title('Señal Original vs PAM Instantáneo');
legend('Original', 'PAM Instantáneo');
grid on;

subplot(3,1,2);
stairs(t, m_pcm, 'r', 'LineWidth', 1.5); % Señal PCM
hold on;
plot(t, m_pam_inst, 'k--');
title(['Señal Cuantificada (N = ', num2str(N), ' bits)']);
legend('Señal Cuantificada', 'PAM Instantáneo');
grid on;

subplot(3,1,3);
histogram(error_quant, 10);
title('Error de cuantificacion');
xlabel('Error');
ylabel('Frecuancia');
grid on;