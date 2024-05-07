function [bw_est, fc_cap]=coarse_BW_fc_estimate(rx_sig, Fs_cor, SNR)

for k=1:length(SNR)
% if SNR>15
%     mid_ind=100;
% elseif SNR>10
%     mid_ind=120;
% else
    mid_ind=140;
%end

%% Estimate noise level

noise = (10^(-SNR(k)/20))*randn(size(rx_sig));

noise_PSD = var(fft(noise));

%% Signal PSD with respect to noise level
sig_fft = fft(rx_sig);
sig_PSD_tmp = 10*log10(abs(sig_fft).^2/noise_PSD);

sig_PSD = 10*log10(smooth((abs(sig_fft).^2)/noise_PSD,20));
sig_PSD(1:20) = 0; % to ignore edge effects in the initial samples

max_snr = sig_PSD(find(sig_PSD==max(sig_PSD)));

detect_threshold = max_snr(1)/1.5;

index = find(sig_PSD(1:floor(end/2)) > detect_threshold);

if length(index)<10;
alg_2 = 0
fc_est=0;
bw_est=0;
else 
alg_2 = 1;

index_diff = diff(index);

ser_index = find(index_diff > mid_ind);

lower = find(ser_index < length(index)/2);
upper = find(ser_index > length(index)/2);

if length(lower) ~= 0 && length(upper) == 0
    index(1:ser_index(lower(end))+1)=[];

elseif length(lower) == 0 && length(upper) ~= 0
    index(ser_index(upper(1))-1:end)=[];

elseif length(lower) ~= 0 && length(upper) ~= 0
    index([1:ser_index(lower(end))+1, ser_index(upper(1))-1:end])=[];
end

fc_cap = ((index(1) + index(end))/2) * Fs_cor/length(sig_PSD);

bw_est = (index(end) - index(1)) * Fs_cor/length(sig_PSD);
%% interpolation 

% low_ind=index(1);
% up_ind=index(end);
% 
% cut_FFT_signal=sig_PSD(low_ind+1:up_ind);
% 
%  f1=(0:length(cut_FFT_signal)-1)*(up_ind*Fs_cor/length(sig_PSD)-(low_ind)*Fs_cor/length(sig_PSD))/length(cut_FFT_signal);
%  f2=(0:I_points-1)*(up_ind*Fs_cor/length(sig_PSD)-(low_ind)*Fs_cor/length(sig_PSD))/I_points;
%  %y=interp1(f1,cut_FFT_signal, f2);
%  
% fc_est2=(low_ind)*Fs_cor/length(sig_PSD)+(I_points/2)*(up_ind*Fs_cor/length(sig_PSD)-(low_ind)*Fs_cor/length(sig_PSD))/I_points;
% nfft = length(sig_PSD);
end
end
 

