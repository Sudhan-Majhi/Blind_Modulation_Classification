clc
clear all
close all

K=1024;           % Number of subcarriers  
n_symbol = 25;    % Number of OFDM symbols
Tuse=3.2e-6;      % useful OFDM symbol period
br=1/Tuse;
Bw=1/Tuse*K;      % delta_f=1/Nuse=carrier pacing 
DT=Tuse/K;        % baseband elementary period=1/sampling rate without oversampling 
G=1/4;            % choice of 1/4, 1/8, 1/16, and 1/32
Tcp=G*Tuse;       % guard band duration
Ts=Tuse+Tcp;      % total OFDM symbol period
SNR=20;
L1=4;             % no. of multipath
oversamp = 50;    % Oversampling factor>20
itr = 10;
Nuse=K*oversamp;
Ncp=K*G*oversamp;
Ns=Nuse+Ncp;
deci = 1;     
count_snr = [];
C_avg = [];
epsilon = 0.5 ;     % Normalized frequency offset
timing_off = 30000; % Timing offset 
C_tile_42_value11 = [];
C_tile_42_value22 = [];
C_tile_42_value33 = [];

for k=1:length(SNR)
for j=1:itr
    
%    Mary={'BPSK','QPSK','PI/4QPSK','MSK','16QAM','OQPSK','8PSK'} ;
     Mary ={'BPSK'} ;
        
      for l=1:length(Mary)            
   
h11 = 1/sqrt(2)*1/sqrt(L1)*(randn(1,L1)+sqrt(-1)*randn(1,L1));    
[tx_sig_base_ofdm] = generate_tx_base_fad_OFDM(K,n_symbol,Ns,oversamp,Mary{l},SNR,h11,Nuse,Ncp) ;      

Dev = randi([-32000 32000],1,n_symbol);

length_Rx_signal = [0:length(tx_sig_base_ofdm)-1] ;
C_epsln_temp_1 = exp(1j*2*pi*(epsilon)/K.*length_Rx_signal);  
tx_sig_base_ofdm = tx_sig_base_ofdm.*C_epsln_temp_1  ;

Ns_est = Ns; % Estimated Ns

tx_sig_base_ofdm = tx_sig_base_ofdm(timing_off:end); % tx_sig_base_ofdm((mm-1)*Ns_est+1+Dev(mm):(mm-1+1)*Ns_est+Dev(mm))  ;tx_sig_base_ofdm(Ns_est+1+Dev(1):2*Ns_est+Dev(1)); tx_sig_base_ofdm(2*Ns_est+1+Dev(2):3*Ns_est+Dev(2));tx_sig_base_ofdm(3*Ns_est+1+Dev(2):4*Ns_est+Dev(2)) ];

          C_tile_42_first = [];
          C_tile_42_second = [];
          C_tile_42_third = [];
          
          for mm = 1:n_symbol-2
    
tx_sig_base_ofdm1 = tx_sig_base_ofdm(mm*Ns_est+1+Dev(mm):(mm+1)*Ns_est+Dev(mm));

%%  Take FFT, 4th EC, for 16QAM, OQPSK, MSK and rest (BPSK, QPSK, PI/4 QPSK) Classification

Rx_sig_base_ofdm_d  = tx_sig_base_ofdm1 ;
rx_OFDM_d = fft(Rx_sig_base_ofdm_d);
    
L0_d = length(rx_OFDM_d);

C_cap_21_d(l) = sum((abs(rx_OFDM_d)).^2)/L0_d;
C_cap_20_d(l) = sum((rx_OFDM_d).^2)/L0_d;
    
% 4th order cumulant 

C_cap_40_d(l) = (sum((rx_OFDM_d).^4)/L0_d)-3*C_cap_20_d(l).^2;
C_cap_41_d(l) = (sum((rx_OFDM_d).^3.*conj(rx_OFDM_d))/L0_d)-3*C_cap_20_d(l)*C_cap_21_d(l);
C_cap_42_d(l) = ((sum((abs(rx_OFDM_d)).^4)/L0_d)-(abs(C_cap_20_d(l))).^2-2*(C_cap_21_d(l)).^2);

% normalized 4th order cumulant

noise1_d = (10^(-SNR/20))*randn(size(rx_OFDM_d))+1i*(10^(-SNR/20))*randn(size(rx_OFDM_d)); 
C_cap_21_noise_d = var(noise1_d);
Rep_C_cap_21_d(l) = C_cap_21_d(l)-C_cap_21_noise_d;

C_tile_42_d(l) = abs(C_cap_42_d(l)/Rep_C_cap_21_d(l).^2) ;

C_tile_42_first = [C_tile_42_first  C_tile_42_d(l)];

%%  Take FFT.^2 and than apply 4th EC, for BPSK, and rest (QPSK, PI/4 QPSK) Classification

Rx_sig_base_ofdm  = tx_sig_base_ofdm1  ;
rx_OFDM = fft(Rx_sig_base_ofdm.^2);
        
rx_fft = rx_OFDM ;
L0 = length(rx_fft);

C_cap_21_fft(l) = sum((abs(rx_fft)).^2)/L0;
C_cap_20_fft(l) = sum((rx_fft).^2)/L0;
    
% 4th order cumulant 

C_cap_40_fft(l) = (sum((rx_fft).^4)/L0)-3*C_cap_20_fft(l).^2;
C_cap_41_fft(l) = (sum((rx_fft).^3.*conj(rx_fft))/L0)-3*C_cap_20_fft(l)*C_cap_21_fft(l);
C_cap_42_fft(l) = ((sum((abs(rx_fft)).^4)/L0)-(abs(C_cap_20_fft(l))).^2-2*(C_cap_21_fft(l)).^2);

% normalized 4th order cumulant

noise1_fft = (10^(-SNR/20))*randn(size(rx_fft))+1i*(10^(-SNR/20))*randn(size(rx_fft)); 
C_cap_21_noise_fft = var(noise1_fft);
Rep_C_cap_21_fft(l) = C_cap_21_fft(l)-C_cap_21_noise_fft;

C_tile_42_fft(l) = abs(C_cap_42_fft(l)/Rep_C_cap_21_fft(l).^2)  ;

C_tile_42_second = [C_tile_42_second  C_tile_42_fft(l)];

%%  Take FFT, than apply DWT(().^2), db45, 4th EC, for QPSK and PI/4 QPSK Classification

Rx_sig_base_db  = fft(tx_sig_base_ofdm1) ;

[corr_seq, lags] = crosscorr(Rx_sig_base_db.^2, conj(Rx_sig_base_db).^2);
% stem(lags, abs(corr_seq)) ;  
% Find strongest magnitude
[Max_Mag_3, in_3] = max(abs(corr_seq));

C_tile_42_third = [C_tile_42_third  Max_Mag_3];

          end
                   
 C_tile_42_value1(l) = sum(C_tile_42_first)/(n_symbol-2);
 C_tile_42_value11 = [C_tile_42_value11 C_tile_42_value1];
 
 C_tile_42_value2(l) = sum(C_tile_42_second)/(n_symbol-2);
 C_tile_42_value22 = [C_tile_42_value22 C_tile_42_value2];

 C_tile_42_value3(l) = sum(C_tile_42_third);
 C_tile_42_value33 = [C_tile_42_value33 C_tile_42_value3];

 %% Algorithm
 
%  [mod_class] = mod_class_algorithm(oversamp,C_tile_42_d,C_tile_42_fft,C_tile_42_dwt_db) ;

            if (C_tile_42_value1 > 88)   
              mod_class='OQPSK';                    
           elseif  (C_tile_42_value1 < 57) 
                mod_class='MSK';                           
           elseif   (81 < C_tile_42_value1 ) && (C_tile_42_value1 <= 88) 
                mod_class='16QAM';                     
           elseif   (57 <= C_tile_42_value1 ) && (C_tile_42_value1 <= 81) &&  (C_tile_42_value2 > 72) 
                mod_class='BPSK';               
           elseif   (57 <= C_tile_42_value1 ) && (C_tile_42_value1 <= 81) &&  (C_tile_42_value2 <= 72) && (C_tile_42_value3 <= 111 )
                mod_class='QPSK';           
           else (57 <= C_tile_42_value1 ) && (C_tile_42_value1 <= 81) &&  (C_tile_42_value2 <= 72) && (C_tile_42_value3 > 111 );
               mod_class='PI/4QPSK';
           end
            
            mod_class          
            count(k,j)=strcmp(mod_class,Mary(:,l))  ;
      end
           
end
snr_count_avg(k)=sum(count(k,:))/itr;
count_snr=[count_snr snr_count_avg(k)];   
 end