clear all; close all; clc

%% load all the data
adcData_breathing_heart_r1 = load('slow_walking_random.mat').adcData_breathing_heart;
adcData_breathing_heart_r2 = load('slow_walking_random2.mat').adcData_breathing_heart;
adcData_breathing_heart_t1 = load('slow_walking_towards.mat').adcData_breathing_heart;
adcData_breathing_heart_t2 = load('slow_walking_towards2.mat').adcData_breathing_heart;
adcData_breathing_heart_t3 = load('slow_walking_towards3.mat').adcData_breathing_heart;
adcData_breathing_heart_t4 = load('slow_walking_towards4.mat').adcData_breathing_heart;

%% view all the data to see what is going on

plot(real(adcData_breathing_heart_r1(:,:,1)))

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

data_single_RX = adcData_breathing_heart_r1(:,:,1);
f = zeros(size(data_single_RX,1), N_samples);

for ii=1:size(data_single_RX,1)
    f(ii,:) = fft(data_single_RX(ii,:),N_samples);
end
% figure()
imagesc(abs(f))