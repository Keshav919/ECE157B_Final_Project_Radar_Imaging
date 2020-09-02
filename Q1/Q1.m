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
    
    % rangeFFT
    RangeFFT = fft(adc_sampled,N_sample,2);
    ant = 1;
    f_max = 1/Ts; %Maximum frequency we can estimate
    tau_max = f_max/k;
    tau_range = 0:tau_resolution:tau_max-tau_resolution;
    distance_range = c*tau_range/2;
    figure; 
    plot(abs(mean(diff(RangeFFT(:,:,ant),1),1)))
    %[~, peak] = max(abs(mean(diff(RangeFFT(:,:,ant),1),1)));
    %fprintf("Dataset "+i+". The peak is at index "+peak+".\n")
    
    % angle range
    theta_range = asind(lambda/2/sRx);%From the lecture
    theta_resolution = asind(lambda/nRx/sRx);%From the lecture

    Nfft = 64; % Number of FFT points for AngleFFT
    u = (-Nfft/2:Nfft/2-1)/Nfft*min([1,2*sin(theta_range)]);% x-axis in world of sin(theta)*(sep/lambda)
    angle_vals = asind(u/sRx*lambda); % theta

    AngleRangeFFT = fftshift(fft(RangeFFT,Nfft,3),3);
    figure; 
    absdiff_anglerange = squeeze(abs(mean(diff(AngleRangeFFT,1,1),1)))';
    surface(distance_range,angle_vals,absdiff_anglerange);
    xlim([0 10])
    xlabel("Distance")
    ylabel("Angle")
    title(["AngleRangeFFT" "_" num2str(i)])
    
    [angle_idx,d_idx] = find(absdiff_anglerange == max(absdiff_anglerange,[],'all'));
    fprintf("In Dataset "+i+", the person is at the location with " + distance_range(d_idx) + "m distance and "+ angle_vals(angle_idx) + " degree angle.\n")
    
    bhFFT = fftshift(fft(unwrap(angle(AngleRangeFFT(1:end-mod(size(RangeFFT,1),2),d_idx,angle_idx)))));
    f = 1/T_frame*(-length(bhFFT)/2:length(bhFFT)/2-1)/length(bhFFT);
    figure
    stem(f,10*log10(abs(bhFFT)))
    xlabel("Frequency (Hz)")
    title(["Phase FFT to get Breathing and HeartBeat" "_" num2str(i)])
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

    fprintf("The breathing rate is "+f(b_idx)+"Hz. The heart rate is "+f(h_idx)+"Hz.\n\n")
end