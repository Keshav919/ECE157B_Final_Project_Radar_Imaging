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

load("two_people.mat");
adc_sampled = adcData_breathing_heart;
%find peak from the rangeFFT
RangeFFT = fft(adc_sampled,N_sample,2);
f_max = 1/Ts; %Maximum frequency we can estimate
tau_max = f_max/k;
tau_range = 0:tau_resolution:tau_max-tau_resolution;
distance_range = c*tau_range/2;
ant = 4;
figure;
plot(abs(mean(diff(RangeFFT(:,:,ant),1),1)))

theta_range = asind(lambda/2/sRx);%From the lecture
theta_resolution = asind(lambda/nRx/sRx);%From the lecture

Nfft = 64; % Number of FFT points for AngleFFT
u = (-Nfft/2:Nfft/2-1)/Nfft*min([1,2*sin(theta_range)]);%% x-axis in world of sin(theta)*(sep/lambda)
angle_vals = asind(u/sRx*lambda);%%x--axis for theta

AngleRangeFFT = fftshift(fft(RangeFFT,Nfft,3),3);
figure; 
%[X,Y] = meshgrid(angle_vals,t*c/2);
absdiff_anglerange = squeeze(abs(mean(diff(AngleRangeFFT,1,1),1)))';
surface(distance_range,angle_vals,absdiff_anglerange);
xlim([0 10])
xlabel("Distance")
ylabel("Angle")

% person 1 location
amp = absdiff_anglerange;
amp(amp<100) = 0;
[x1,y1] = find(imregionalmax(amp)==1);

bhFFT = fftshift(fft(unwrap(angle(AngleRangeFFT(1:end-mod(size(RangeFFT,1),2),x1,y1))))); % angle based analyzing
f = 1/T_frame*(-length(bhFFT)/2:length(bhFFT)/2-1)/length(bhFFT);
figure
stem(f,10*log10(abs(bhFFT)))
xlim([-2,2])

% breathing peak finding
bfiltered = abs(bhFFT);
bfiltered((f<0.1) | (f>0.5)) = 0;
[~, b_unfilt_idx] = findpeaks(bfiltered);
b_filtered_idx = b_unfilt_idx((b_unfilt_idx>max(find(f<0.21))) & (b_unfilt_idx<min(find(f>0.42))));
[~, b_filt_idx] = max(bfiltered(b_filtered_idx));
b_idx = b_filtered_idx(b_filt_idx);
% heart peak finding
hfiltered = abs(bhFFT);
hfiltered((f<0.9) | (f>1.6)) = 0;
[~, h_unfilt_idx] = findpeaks(hfiltered);
h_filtered_idx = h_unfilt_idx((h_unfilt_idx>max(find(f<1.1))) & (h_unfilt_idx<min(find(f>1.5))));
[~, h_filt_idx] = max(hfiltered(h_filtered_idx));
h_idx = h_filtered_idx(h_filt_idx);


fprintf("The breathing rate is "+f(b_idx)+"Hz. The heart rate is "+f(h_idx)+"Hz.\n")


% person 2 location and breathing and heart rate
amp = absdiff_anglerange;
amp(:,y1-5:y1+5) = 0;
amp(amp<40) = 0;
[x2,y2] = find(imregionalmax(amp)==1);

bhFFT = fftshift(fft(unwrap(angle(AngleRangeFFT(1:end-mod(size(RangeFFT,1),2),x2,y2))))); % angle based analyzing
f = 1/T_frame*(-length(bhFFT)/2:length(bhFFT)/2-1)/length(bhFFT);
figure
stem(f,10*log10(abs(bhFFT)))
xlim([-2,2])

% breathing peak finding
bfiltered = abs(bhFFT);
bfiltered((f<0.1) | (f>0.5)) = 0;
[~, b_unfilt_idx] = findpeaks(bfiltered);
b_filtered_idx = b_unfilt_idx((b_unfilt_idx>max(find(f<0.21))) & (b_unfilt_idx<min(find(f>0.42))));
[~, b_filt_idx] = max(bfiltered(b_filtered_idx));
b_idx = b_filtered_idx(b_filt_idx);
% heart peak finding
hfiltered = abs(bhFFT);
hfiltered((f<0.9) | (f>1.6)) = 0;
[~, h_unfilt_idx] = findpeaks(hfiltered);
h_filtered_idx = h_unfilt_idx((h_unfilt_idx>max(find(f<1.1))) & (h_unfilt_idx<min(find(f>1.5))));
[~, h_filt_idx] = max(hfiltered(h_filtered_idx));
h_idx = h_filtered_idx(h_filt_idx);

fprintf("The breathing rate is "+f(b_idx)+"Hz. The heart rate is "+f(h_idx)+"Hz.\n")