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

% plot(real(adcData_breathing_heart_r1(:,:,1)))

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

adc_sampled = adcData_breathing_heart_t1(100:end-10,:,:);
RangeFFT = fft(adc_sampled(:,:,ant),N_sample,2);
N_frame = size(adc_sampled,1);
t_frame = (1:size(adc_sampled,1))*T_frame;

f_max = 1/Ts; %Maximum frequency we can estimate
tau_max = f_max/k;
tau_range = 0:tau_resolution:tau_max-tau_resolution;
distance_range = c*tau_range/2;

% figure
% imagesc(distance_range, t_frame, abs(RangeFFT))

% Plot differentiated peaks
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

%% Doppler FFT (range limitation hard to estimate)
%stft(RangeFFT,1/T_frame,'Window',kaiser(256,5),'OverlapLength',220,'FFTLength',512);
% N_frame = N_frame/8;
% 
% DopplerFFT = fft(RangeFFT(1:N_frame,:),N_frame,1);
% % max_index=find(RangeFFT(1,:)==max(RangeFFT(1,:)));
% v_max = (c/(2*f_start))/T_frame;
% v_res = v_max/(N_frame);
% y_axis = v_res:v_res:v_max;
% figure; imagesc(distance_range,y_axis,abs(DopplerFFT));
% xlabel('Range')
% ylabel('doppler')

%v_init = (distance_range(range_idx(500))- distance_range(range_idx(400)))/100/T_frame;

%% 1D Kalmen filter - ref: http://studentdavestutorials.weebly.com/kalman-filter-with-matlab-code.html

dt = T_frame; % sample time
P_loc_meas = distance_range(range_idx); % person path
v_measure = diff(P_loc_meas)/dt;
v_measure = movmean(v_measure,200);
P_loc_meas = P_loc_meas(101:end-101);

% Define update equations (Coefficent matrices)
A = [1 dt; 0 1] ; % state transition matrix
B = [dt^2/2; dt]; % input control matrix
C = [1 0; 0 1]; % measurement matrix, only measuring position

% Define main variables
u = 0; % define acceleration magnitude, constant velocity
X = [P_loc_meas(1); -0.06]; % initized state - [position; velocity]
X_estimate = X;  % x_estimate of initial location estimation
personAccel_noise_mag = 1e-3; % process noise
Measure_noise_x = 0.05; % location measurement noise
Measure_noise_v = 0.1; % velocity measurement noise
R = [Measure_noise_x^2 0; 0 Measure_noise_v^2];% Ez convert the measurement noise (stdv) into covariance matrix
Ex = personAccel_noise_mag^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial person position variance (covariance matrix)

%initize estimation variables
P_loc_estimate = []; % person position estimate
vel_estimate = []; % person velocity estimate
P_mag_estimate = [];

for t = 1:length(P_loc_meas)
    % Predict next state of the person with the last state and predicted motion.
    X_estimate = A * X_estimate + B * u;
    % predict next covariance
    P = A * P * A'+ Ex;
    % Kalman Gain
    K = P*C'*inv(C*P*C'+R);
    % Update the state estimate.
    X_estimate = X_estimate + K * ([P_loc_meas(t); v_measure(t)] - C * X_estimate);
    % Update covariance estimation.
    P =  (eye(2)-K*C)*P;
    % Store for plotting
    P_loc_estimate = [P_loc_estimate; X_estimate(1)];
    vel_estimate = [vel_estimate; X_estimate(2)];
    P_mag_estimate = [P_mag_estimate; P(1)];
end
hold on;
plot(P_loc_estimate,t_frame(101:end-101),'r')
legend('estimated')

%% Breathing

% [distance, est] = meshgrid(distance_range,P_loc_estimate);
% est_idx = sum(distance_range<P_loc_estimate,2); % remapping the index is not quite useful 

% Compensate the phase change introduced by distance
angle_d_compensated = zeros(1,2*floor((N_frame-101)/2)-100);
for i = 101:2*floor((N_frame-101)/2)
    angle_d_compensated(i-100) = angle(RangeFFT(i,range_idx(i))) + 2*pi*f_start/c*(P_loc_estimate(i-100));
end

figure
plot(unwrap(angle_d_compensated))

bhFFT = fftshift(fft(unwrap(angle_d_compensated)));
f = 1/T_frame*(-length(bhFFT)/2:length(bhFFT)/2-1)/length(bhFFT);

figure
stem(f,10*log10(abs(bhFFT)))
xlim([-2,2])

% Breathing peak finding
bfiltered = abs(bhFFT);
bfiltered((f<0.1) | (f>0.5)) = 0;
[~, b_unfilt_idx] = findpeaks(bfiltered);
b_filtered_idx = b_unfilt_idx((b_unfilt_idx>max(find(f<0.23))) & (b_unfilt_idx<min(find(f>0.42))));
[~, b_filt_idx] = max(bfiltered(b_filtered_idx));
b_idx = b_filtered_idx(b_filt_idx);
fprintf("The breathing rate is "+f(b_idx)+"Hz.\n")