clear all; close all; clc

%% load all the data
adcData_breathing_heart = load('adcData_Exercising.mat').adcData_breathing_heart;
adcData_breathing_heart_s1 = load('adcData_standing_exercising1.mat').adcData_breathing_heart;
adcData_breathing_heart_s2 = load('adcData_standing_exercising2.mat').adcData_breathing_heart;
adcData_breathing_heart_s3 = load('adcData_standing_exercising3.mat').adcData_breathing_heart;

%% view all the data to see what is going on

plot(real(adcData_breathing_heart(:,:,1)))

%% constants

f_start = 77e9;
k = 79;
t_chirp = 10e-3;
fs = 1/t_chirp;

N_samples = 256;

fs_adc = 8000e3;
ts_adc = 1/fs_adc;
M = 3000;

Nrx = 4;

c=3e8;


%% Range FFT

data_single_RX = adcData_breathing_heart(:,:,1);
f = zeros(size(data_single_RX,1), N_samples);

for ii=1:size(data_single_RX,1)
    f(ii,:) = fft(data_single_RX(ii,:),N_samples);
end
% figure()
imagesc(abs(f))

%% Analyse Phase
f_mean = mean(abs(f));
[~,max_indx] = max(abs(f_mean(:)));

person = f(:,max_indx);
phase = unwrap(angle(person));
plot(abs(phase))
%%
plot(abs(f_mean))
%%

fft_len = length(phase)*3;
freq_axis = (-fft_len/2:fft_len/2-1)*fs/size(data_single_RX,1);

fbhh = fftshift(fft(phase, fft_len));


plot(freq_axis,abs(fbhh));
xlabel("Frequency (Hz)")
xlim([-2, 2])
title("Frequency analysis of Human")
