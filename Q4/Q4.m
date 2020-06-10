clear all; close all; clc

%% load all the data
load('slow_walking_random.mat')
adcData_breathing_heart_r1 = adcData_breathing_heart;
load('slow_walking_random2.mat')
adcData_breathing_heart_r2 = adcData_breathing_heart;
load('slow_walking_towards.mat')
adcData_breathing_heart_t1 = adcData_breathing_heart;
load('slow_walking_towards2.mat')
adcData_breathing_heart_t2 = adcData_breathing_heart;
load('slow_walking_towards3.mat')
adcData_breathing_heart_t3 = adcData_breathing_heart;
load('slow_walking_towards4.mat')
adcData_breathing_heart_t4 = adcData_breathing_heart;

%% view all the data to see what is going on

plot(real(adcData_breathing_heart_r1(:,:,1)))

%% constants
f_start = 77e9;

k = 79e12;
c = 3e8;

N_sample = 256;
fs = 8e6;
Ts = 1/fs;
T = N_sample*Ts;

T_frame = 10e-3;
M = 3000;

t_sample = (0:N_sample-1)*Ts;

B = k*N_sample*Ts;
tau_resolution = 1/B;

lambda = c/f_start;
nRx = 4; % number of antennas
sRx = lambda/2;% separation between receivers

ant = 1;

%% Range FFT

adc_sampled = adcData_breathing_heart_r2;
RangeFFT = fft(adc_sampled(:,:,ant),N_sample,2);
N_frame = size(adc_sampled,1);
t_frame = (1:size(adc_sampled,1))*T_frame;

f_max = 1/Ts; %Maximum frequency we can estimate
tau_max = f_max/k;
tau_range = 0:tau_resolution:tau_max-tau_resolution;
distance_range = c*tau_range/2;

figure
imagesc(distance_range, t_frame, abs(RangeFFT))

figure
imagesc(distance_range, t_frame, abs(diff(RangeFFT)))

[~, range_idx] = max(abs(diff(RangeFFT)),[],2);
range_idx = [range_idx(1); range_idx];

% %% Angle FFT
% 
% theta_range = asind(lambda/2/sRx);%From the lecture
% theta_resolution = asind(lambda/nRx/sRx);%From the lecture
% 
% Nfft = 64; % Number of FFT points for AngleFFT
% u = (-Nfft/2:Nfft/2-1)/Nfft*min([1,2*sin(theta_range)]);%% x-axis in world of sin(theta)*(sep/lambda)
% angle_vals = asind(u/sRx*lambda);%%x--axis for theta
% 
% RangeFFT_rmrange = zeros(size(adc_sampled,1), size(adc_sampled,3));
% for i = 1:size(adc_sampled,1)
%     RangeFFT_rmrange = squeeze(RangeFFT(i,range_idx,:));
% end
% 
% AngleRangeFFT = fftshift(fft(RangeFFT_rmrange,Nfft,2),2);
% figure; 
% surface(angle_vals,t_frame,abs(AngleRangeFFT),'EdgeColor','none');
% %xlim([0 10])
% xlabel("Time")
% ylabel("Angle")
% 
% [~, angle_idx] = max(abs(AngleRangeFFT),[],2);

%% Doppler FFT

DopplerFFT = fft(RangeFFT,N_frame,1);
max_index=find(RangeFFT(1,:)==max(RangeFFT(1,:)));
v_max = (c/(2*f_start))*(f_breathing_samp);
v_res = v_max/(N_chirps);
y_axis = v_res:v_res:v_max;
figure; imagesc(distance_range,y_axis,abs(DopplerFFT));
xlabel('Range')
ylabel('doppler')
