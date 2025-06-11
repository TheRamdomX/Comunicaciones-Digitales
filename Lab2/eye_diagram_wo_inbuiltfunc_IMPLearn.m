% Code made by Manu Prasad (https://www.youtube.com/@IMPLearn)
% Program for plotting eye diagram of pulse shaped signal using raised
% cosine pulse without using inbuilt functions
clc; 
clear; 
close all; 
N  = 300; % number of symbols 
data = randi([0,1],1,N); 
am = 2.*data -1;% converting to polar +1 & -1 
fs = 10; % sampling rate
t = -fs:1/fs:fs;   

% defining the sinc filter 
sincNum = sin(pi*t);    % numerator of the sinc function 
sincDen = (pi*t);       % denominator of the sinc function 
sincDenZero = abs(sincDen) < 10^-10; % Finding index of values closer to zero
sincOp = sincNum./sincDen; 
sincOp(sincDenZero) = 1; % sin(pix)/(pix) =1 for x =0   

% raised cosine filter for alpha = 0.5 
alpha = 0.5; 
cosNum = cos(alpha*pi*t); 
cosDen = (1-(2*alpha*t).^2); 
cosDenZero = abs(cosDen)<10^-10;% Finding index of values closer to zero
cosOp = cosNum./cosDen; 
cosOp(cosDenZero) = pi/4;   
gt_alpha5 = sincOp.*cosOp;   
% Plotting 
figure; 
plot(t,gt_alpha5); 
title('Raised Cosine pulse for alpha = 0.5');   

% raised cosine filter for alpha = 1 
alpha = 1; 
cosNum = cos(alpha*pi*t); 
cosDen = (1-(2*alpha*t).^2); 
cosDenZero = find(abs(cosDen)<10^-10); 
cosOp = cosNum./cosDen; 
cosOp(cosDenZero) = pi/4;   
gt_alpha1 = sincOp.*cosOp;   
% Plotting
figure; 
plot(t,gt_alpha1); 
title('Raised Cosine pulse for alpha = 1');   

% upsampling the transmit sequence Method 1
% amUpSampled = [am;zeros(fs-1,length(am))]; 
% amU = amUpSampled(:).'; 

% upsampling the transmit sequence Method 2
am_up = upsample(am,fs);   

% filtered sequence Method 1
% st_alpha5 = conv(amU,gt_alpha5); 
% st_alpha1 = conv(amU,gt_alpha1);   

% filtered sequence Method 2
st_alpha5 = conv(am_up,gt_alpha5); 
st_alpha1 = conv(am_up,gt_alpha1);   

% taking only the first N*fs samples  
st_alpha5_rs = st_alpha5(1:(N*fs)); 
st_alpha1_rs = st_alpha1(1:(N*fs));  

% Reshaping the signal to [fs*2 by N/2]
st_alpha5_reshape = reshape(st_alpha5_rs,fs*2,N/2); 
st_alpha1_reshape = reshape(st_alpha1_rs,fs*2,N/2);   

% Plotting eye diagram for alpha =0.5 
t1 = 0:1/fs:(((2*fs)-1)/fs); 
figure; 
plot(t1,real(st_alpha5_reshape),'b');    
title('eye diagram with alpha=0.5'); 
xlabel('time') 
ylabel('amplitude')  
axis([0 2 -1.5 1.5]) 
grid on   
% Plotting eye diagram for alpha =1 
figure; 
plot(t1,real(st_alpha1_reshape),'b');   
title('eye diagram with alpha=1')  
xlabel('time')
ylabel('amplitude')  
axis([0 2 -1.5 1.5 ]) 
grid on