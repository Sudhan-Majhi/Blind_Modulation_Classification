 clear all
 close all
 clc

% tic             % Start of stopwatch timer (toc at bottom)

L = 1000;        % Number of symbols
R = 0.2;          % Roll-off factor for root raised cosine (RRC) filter
delay = 3;        % Delay of RRC filter
             
SNR=0:5:20;       % Receiver average SNR
SNR= 15;
nTx = 2;          % number of tx antennas  
nRx = 2;          % number of rx antennas 
fs = 1e6;         % Symbol rate [Hz]
Fs = 45e6;        % Sampling frequency [Hz]
sim_Fs = fs*820;
itr = 10;          % no. of iterations
L1 = 4;           % no. of multipaths
oversamp = 10;    % carrier period to elementary period ratio >5  

deci =2;
Fs_cor = oversamp*fs/deci;     %this should be eight times the carrier frequency for our algo to estimate without ambiguity
theta=0;
fe=0;
the=0;
count_snr=[];
%Mary={'QPSK','OQPSK','PI/4QPSK','MSK','8PSK','16QAM'};
Mary ={'MSK'};
    
for k=1:length(SNR) 
  M=[];  
 for j=1:itr % itr=5   

  for l=1:length(Mary) 
%% Transmitter Modulation schemes 
%  Multipath channel coefficients 

h11 = 1/sqrt(2)*1/sqrt(L1)*(randn(1,L1)+sqrt(-1)*randn(1,L1));
h12 = 1/sqrt(2)*1/sqrt(L1)*(randn(1,L1)+sqrt(-1)*randn(1,L1));
h21 = 1/sqrt(2)*1/sqrt(L1)*(randn(1,L1)+sqrt(-1)*randn(1,L1));
h22 = 1/sqrt(2)*1/sqrt(L1)*(randn(1,L1)+sqrt(-1)*randn(1,L1));
% RRC pulse
rc = rcosine(1,oversamp,'sqrt',R);  

 [tx_sig_1,tx_sig_2] = generate_sg_fad_baseband_MIMO(L,oversamp,deci,Fs_cor,Mary{l},theta,SNR,h11,h12,h21,h22,delay,rc,nTx);


%% BW estimation @ antenna 1
          
      [bw_est_1] = coarse_BW_fc_estimate(tx_sig_1, Fs_cor, SNR);      
%        bw_est_1=2385000;
      low_fs = bw_est_1*.30;    
      up_fs = bw_est_1*.80;      
      low_ind1 = round(low_fs*length(tx_sig_1)/Fs_cor);
      up_ind1 = round(up_fs*length(tx_sig_1)/Fs_cor);
  
      %% BW estimation @ antenna 2
                 
       [bw_est_2] = coarse_BW_fc_estimate(tx_sig_2, Fs_cor, SNR);      
%        bw_est_2=2385000;
      low_fs = bw_est_2*.30;    
      up_fs = bw_est_2*.80;      
      low_ind2 = round(low_fs*length(tx_sig_2)/Fs_cor);
      up_ind2 = round(up_fs*length(tx_sig_2)/Fs_cor);
      
%% 2nd order cyclic cumulant at baseband level at antenna 1      

     rx_bb_1 = abs(tx_sig_1).^2;   %  ||^2             
    rx_bb_1 = rx_bb_1/norm(rx_bb_1);  
     f = (0:length(rx_bb_1)/2-1)*Fs_cor/length(rx_bb_1);
     fft_sig_1 = abs(fft(rx_bb_1));      %     DFT
     fft_sig_1 = fft_sig_1(1:length(rx_bb_1)/2);
     
%% 2nd order cyclic cumulant at baseband level at antenna 2     

    rx_bb_2 = abs(tx_sig_2).^2;   %  ||^2              
    rx_bb_2 = rx_bb_2/norm(rx_bb_2);
    fft_sig_2 = abs(fft(rx_bb_2.^2));      %     DFT
    fft_sig_2 = fft_sig_2(1:length(rx_bb_2)/2);
    
 %% 2nd order Cyclic peak QPSK, PI/4QPSK, MSK analysis at ws and 2ws
  
%   figure(2)
    hold on
    %plot(f,fft_sig_1,'bd-','Linewidth',1) 

    fft_sig_1(1:low_ind1)=0;
%     fft_sig_1(4*up_ind1+1:end)=0;      
 
%   figure(6)
    plot(f,fft_sig_1(1:length(fft_sig_1)))
    [maxQ1, index1] = max(fft_sig_1);
    f_sym_peak_Q_1 = f(index1);
 
%   figure(7)
%   plot(f,fft_sig_1(1:length(rx_bb_2)/2))        
     fft_sig_2(1:low_ind2)=0;
%      fft_sig_2(4*up_ind2+1:end)=0;      
    
%   figure(8)
%   plot(f,fft_sig_2(1:length(fft_sig_2)/2)) 
    [maxQ2, index2] = max(fft_sig_2);
    f_sym_peak_Q_2 = f(index2);


%% 2nd order correlation 
% tx_sig_1 = tx_sig_1/norm(tx_sig_1);
% tx_sig_2 = tx_sig_2/norm(tx_sig_2);
  [corr_seq, lags]= crosscorr(tx_sig_1, conj(tx_sig_2));
%   dfittool(abs(corr_seq))
%    figure(13)
% hold on
     stem(lags, abs(corr_seq),'k');  
%   Find strongest magnitude
  [Max_Mag_1, in_1] = max(abs(corr_seq));
Max_Mag_1;
%% 4th order correlation 
 
  [corr_seq, lags] = crosscorr(tx_sig_1.^2, conj(tx_sig_2).^2);
%   figure(14)
%    hold on
     stem(lags, abs(corr_seq),'R-s') ;  
%   Find strongest magnitude
    [Max_Mag_2, in_2] = max(abs(corr_seq));
Max_Mag_2;


%% 4th order correlation  for classifying 8PSK and PI/4QPSK, we convert PI/4QPSK to QPSK to get the feature (peak for PI/4QPSK)
  
 tx_sig_11 = tx_sig_1(1:oversamp:end);
 
 tx_sig_22 = tx_sig_2(1:oversamp:end);
 
  tx_sig_11 = tx_sig_11(1:2:end);
  tx_sig_11 = repmat(tx_sig_11,1,oversamp);
  
  tx_sig_22 = tx_sig_22(1:2:end);
  tx_sig_22 = repmat(tx_sig_22,1,oversamp);
  
  [corr_seq, lags] = crosscorr(tx_sig_11.^2, conj(tx_sig_22).^2);
%   figure(15)
%   hold on
 stem(lags, abs(corr_seq)) ;  
%   Find strongest magnitude
    [Max_Mag_3, in_3] = max(abs(corr_seq));
Max_Mag_3;

%% 4th elementary cumulant for classifying QPSK and 16QAM

L0 = length(tx_sig_1);    % at antenna 1
C_cap_211=sum((abs(tx_sig_1)).^2)/L0;
C_cap_201=sum((tx_sig_1).^2)/L0;


%4th order cumulant 
C_cap_401=(sum((tx_sig_1).^4)/L0)-3*C_cap_201.^2;
C_cap_411=(sum((tx_sig_1).^3.*conj(tx_sig_1))/L0)-3*C_cap_201*C_cap_211;
C_cap_421 =((sum((abs(tx_sig_1)).^4)/L0)-(abs(C_cap_201)).^2-2*(C_cap_211).^2);

C_cap_212=sum((abs(tx_sig_2)).^2)/L0;   % at antenna 2
C_cap_202=sum((tx_sig_2).^2)/L0;


%4th order cumulant 
C_cap_402=(sum((tx_sig_2).^4)/L0)-3*C_cap_202.^2;
C_cap_412=(sum((tx_sig_2).^3.*conj(tx_sig_2))/L0)-3*C_cap_202*C_cap_212;
C_cap_422=((sum((abs(tx_sig_2)).^4)/L0)-(abs(C_cap_202)).^2-2*(C_cap_212).^2);

C_cap_42_avg = (C_cap_422+C_cap_421)/2;



  end
 end
end