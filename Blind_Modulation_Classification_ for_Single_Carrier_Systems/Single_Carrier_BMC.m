clear all
close all
clc
fc = 5e6;                     % Carrier frequency [Hz]
fc_search = 5;                % fc belongs within +-5% of coarse fc
Fs = 45e6;                    % Sampling frequency [Hz]
fs = 122e3;                   % Symbol rate [Hz]
sim_Fs = fs*820;
L = 200;                      % Number of symbols
R = .5;                       % Roll-off factor for root raised cosine (RRC) filter
delay=3;                      % Delay of RRC filter
SNR=-5:5:15;                  % Receiver average SNR
SNR=20;
itr = 1;                      % Number of iterations
T = 1/Fs;                     % Sample time
t =(0:L-1)*T;                 % Time vector
oversamp = 820;             % Oversampling factor
deci = (sim_Fs - rem(sim_Fs,Fs))/Fs;    % Non-integer
Fs_cor = oversamp*fs/deci;     %this should be eight times the carrier frequency for our algo to estimate without ambiguity

rc = rcosine(1,oversamp,'sqrt',R);
theta=0;
fe=0;
the=0;
%  Mary={'BPSK','QPSK','OQPSK','PI/4QPSK','MSK','16QAM'};
Mary ={'QPSK'};


for k=1:length(SNR)
    for j=1:itr
        for l=1:length(Mary)
            
            [tx_sig] = generate_transmit_signal_fading(L,fc,oversamp,deci,Fs_cor,Mary{l}, fe, theta,R,SNR);
            rx_sig = tx_sig ;
            rx_sig = hilbert(rx_sig) ;
            rx_norm = rx_sig/sqrt(mean(abs(rx_sig).^2)) ;
            
            % cumulant in frequency domain
            f=(0:length(rx_norm)/2-1)*Fs_cor/length(rx_norm);
            
            fft_rx = fft(rx_norm.^4)/length(rx_norm);
            
            fft_rx1 = abs(fft_rx(1:length(rx_sig)/2)).^2;
            %  plot(f,fft_rx1)
            noise = (10^(-SNR(k)/20))*randn(size(rx_norm));
            noise_PSD =sqrt( var(fft(noise)));
            fft_n = fft_rx/noise_PSD ;
            
            fft_rn = abs(fft_n(1:length(rx_sig)/2)).^2;
            %Find strongest magnitude
            [fft_max index] = max(fft_rn);
            
            %Look up the frequency that corresponds with the strongest magnitude
            f_cyc = f(index);
            % plot(f,fft_rx1)
            
            %%  downconversion & remove higher freq comp
            
            time_vec = (0:length(rx_norm)-1)/Fs_cor;
            [bw_est, fc_cap]=coarse_BW_fc_estimate(rx_sig, Fs_cor, SNR);
            fc_cap=5e6;
            
            low_fs=bw_est*.30;
            up_fs=bw_est*.80;
            
            low_ind=round(low_fs*length(rx_norm )/Fs_cor);
            up_ind=round(up_fs*length(rx_norm)/Fs_cor);
            
            I_comp_dc = rx_norm.* cos(2*pi*fc_cap*time_vec);  % I component
            nfft=length(I_comp_dc);
            bn= abs(bw_est*nfft/Fs_cor);
            LPF=round(1+abs(bw_est*nfft/Fs_cor)) ;
            fft_I_comp = fft(I_comp_dc);
            fft_I_comp(LPF:end-LPF) = 0;
            I_rx_sig = ifft(fft_I_comp);
            
            Q_comp_dc = -rx_norm.* sin(2*pi*fc_cap*time_vec);
            
            LPF=round(1+abs(bw_est*nfft/Fs_cor)) ;
            fft_Q_comp = fft(Q_comp_dc);
            fft_Q_comp(LPF:end-LPF) = 0;
            Q_rx_sig = ifft(fft_Q_comp);
            
            % I and Q component with complex representation
            IQ_comp = I_rx_sig+ 1i*Q_rx_sig ;
            rx_bb=IQ_comp;
            
            %% 2nd order cumulant at baseband level
            f=(0:length(rx_bb)/2-1)*Fs_cor/length(rx_bb);
            fft_sig_2=abs(fft(abs(rx_bb).^2));      %     DFT
            
            plot(f,fft_sig_2(1:length(rx_bb)/2))
            
            L0=length(rx_bb);
            
            C_cap_21(l)=sum((abs(rx_bb)).^2)/L0;
            C_cap_20(l)=sum((rx_bb).^2)/L0;
            
            
            %4th order cumulant
            C_cap_40(l)=(sum((rx_bb).^4)/L0)-3*C_cap_20(l).^2;
            C_cap_41(l)=(sum((rx_bb).^3.*conj(rx_bb))/L0)-3*C_cap_20(l)*C_cap_21(l);
            C_cap_42(l)=((sum((abs(rx_bb)).^4)/L0)-(abs(C_cap_20(l))).^2-2*(C_cap_21(l)).^2);
            
            
            %normalized 4th order cumulant
            noise1 = (10^(-SNR(k)/20))*randn(size(tx_sig))+1i*(10^(-SNR(k)/20))*randn(size(tx_sig));
            C_cap_21_noise=var(noise1);
            Rep_C_cap_21=C_cap_21(l)-C_cap_21_noise;
            
            C_tile_40(l)=C_cap_40(l)/Rep_C_cap_21;
            C_tile_41(l)=C_cap_41(l)/Rep_C_cap_21;
            C_tile_42(l)=C_cap_42(l)/Rep_C_cap_21;
               
        end
    end %itr
end %snr
toc



