%% FMCW
clear
close all
f_start = 77e9;

k = 79e12;
c = 3e8;

N_sample = 256;
fs = 8e6;
Ts = 1/fs;
T = N_sample*Ts;

T_frame = 10e-3;
M = 3000;

t = (0:N_sample-1)*Ts;

B = k*N_sample*Ts;
tau_resolution = 1/B;


lambda = c/f_start;
nRx = 4; % number of antennas
sRx = lambda/2;% separation between receivers

for i = 1:4
    load("adcData_Q1data"+i+".mat");
    adc_sampled = adcData_breathing_heart;
    %find peak from the rangeFFT
    RangeFFT = fft(adc_sampled,N_sample,2);%??????
    ant = 1;
    figure; 
    plot(abs(mean(diff(RangeFFT(:,:,ant),1),1)))

    [~, peak] = max(abs(mean(diff(RangeFFT(:,:,ant),1),1)));%??????
    fprintf("The peak is at index"+peak+".\n")

    bhFFT = fftshift(fft(unwrap(angle(RangeFFT(:,peak,2)))));
    f = 1/T_frame*(-length(bhFFT)/2:length(bhFFT)/2-1)/length(bhFFT);
    figure
    stem(f,10*log10(abs(bhFFT)))
    xlim([-2,2])
end