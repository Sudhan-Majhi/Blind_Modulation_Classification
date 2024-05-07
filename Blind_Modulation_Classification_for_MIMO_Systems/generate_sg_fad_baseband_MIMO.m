function [tx_sig_1,tx_sig_2] = generate_sg_fad_baseband_MIMO(L,oversamp,deci,Fs_cor,Mary,theta,SNR,h11,h12,h21,h22,delay,rc,nTx)


for k=1:length(SNR)
    
    switch Mary
        
        case 'BPSK'
            
            x_real = 2*randi(1,L)-1;
            x_filtered = rcosflt(x_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered1 = x_filtered(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            sig = reshape(x_filtered1,length(x_filtered1)/nTx,nTx)';
            
            ant_1 = sig(1,:);    % symbol txd from antenna 1
            ant_2 = sig(2,:);    % symbol txd from antenna 2
            
            ant11_conv = conv(ant_1,h11);   % convolution of each symbol with channel h11
            ant11_conv = ant11_conv(1:length(ant_1));
            
            ant12_conv = conv(ant_2,h12);   % convolution of each symbol with channel h12
            ant12_conv = ant12_conv(1:length(ant_2));
            
            ant21_conv = conv(ant_1,h21);   % convolution of each symbol with channel h21
            ant21_conv = ant21_conv(1:length(ant_1));
            
            ant22_conv = conv(ant_2,h22);   % convolution of each symbol with channel h22
            ant22_conv = ant22_conv(1:length(ant_2));
            
            n1 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 1
            n2 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 2
            
            y1 = ant11_conv + ant12_conv + 10^(-SNR(k)/20)*n1 ;
            y2 = ant21_conv + ant22_conv + 10^(-SNR(k)/20)*n2 ;
            
            y1 = y1(1:deci:end);
            y2 = y2(1:deci:end);
            
            t = (0:length(y1)-1)/Fs_cor;
            
            
            tx_sig_1 = y1 ;  % at antenna 1
            
            tx_sig_2 = y2 ;  % at antenna 2
            
            
        case 'QPSK'
            xa=2*randi(1,2*L)-1;
            
            sig = reshape(xa,[],nTx)';
            ant_1 = sig(1,:);    % symbol txd from antenna 1
            
            
            x_evn=ant_1(1:2:end);
            x_odd=ant_1(2:2:end);
            
            x_filtered_r = rcosflt(x_evn,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r1 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            x_filtered_i = rcosflt(x_odd,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i1 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            signal_1 = x_filtered_r1 + 1j*x_filtered_i1;
            
            ant_2 = sig(2,:);    % symbol txd from antenna 2
            
            x_evn=ant_2(1:2:end);
            x_odd=ant_2(2:2:end);
            
            x_filtered_r = rcosflt(x_evn,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r1 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            x_filtered_i = rcosflt(x_odd,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i1 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            signal_2 = x_filtered_r1 + 1j*x_filtered_i1;
            
            % sig = reshape(signal,length(signal)/nTx,nTx)';
            
            % ant_1 = sig(1,:);    % symbol txd from antenna 1
            % ant_2 = sig(2,:);    % symbol txd from antenna 2
            
            ant11_conv = conv(signal_1,h11);   % convolution of each symbol with channel h11
            ant11_conv = ant11_conv(1:length(signal_1));
            
            ant12_conv = conv(signal_2,h12);   % convolution of each symbol with channel h12
            ant12_conv = ant12_conv(1:length(signal_2));
            
            ant21_conv = conv(signal_1,h21);   % convolution of each symbol with channel h21
            ant21_conv = ant21_conv(1:length(signal_1));
            
            ant22_conv = conv(signal_2,h22);   % convolution of each symbol with channel h22
            ant22_conv = ant22_conv(1:length(signal_2));
            
            n1 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 1
            n2 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 2
            
            y1 = ant11_conv + ant12_conv + 10^(-SNR(k)/20)*n1 ;
            y2 = ant21_conv + ant22_conv + 10^(-SNR(k)/20)*n2 ;
            
            %  y1 = ant11_conv ;
            %  y2 = ant21_conv + ant22_conv ;
            
            y1 = y1(1:deci:end);
            y2 = y2(1:deci:end);
            
            t = (0:length(y1)-1)/Fs_cor;
            
            tx_sig_1 = y1;  % at antenna 1
            
            tx_sig_2 = y2;  % at antenna 2
            
            
        case 'OQPSK'
            
            xa=2*randi(1,2*L)-1;
            
            sig = reshape(xa,[],nTx)';
            ant_1 = sig(1,:);    % symbol txd from antenna 1
            
            e_xa=ant_1(1:2:end);
            o_xa=ant_1(2:2:end);
            
            
            e_filtered = rcosflt(e_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            e_filtered1 = e_filtered(oversamp*delay+1:end-oversamp*delay)';
            x_shift_r = [e_filtered1, zeros(1, floor(oversamp/2))];
            
            o_filtered = rcosflt(o_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            o_filtered1 = o_filtered(oversamp*delay+1:end-oversamp*delay)';
            x_shift_i = [zeros(1, floor(oversamp/2)), o_filtered1];
            
            signal_1 = x_shift_r + 1j*x_shift_i;
            
            ant_2 = sig(2,:);    % symbol txd from antenna 1
            
            e_xa=ant_2(1:2:end);
            o_xa=ant_2(2:2:end);
            
            
            e_filtered = rcosflt(e_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            e_filtered1 = e_filtered(oversamp*delay+1:end-oversamp*delay)';
            x_shift_r = [e_filtered1, zeros(1, floor(oversamp/2))];
            
            o_filtered = rcosflt(o_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            o_filtered1 = o_filtered(oversamp*delay+1:end-oversamp*delay)';
            x_shift_i = [zeros(1, floor(oversamp/2)), o_filtered1];
            
            signal_2 = x_shift_r + 1j*x_shift_i;
            
            % sig = reshape(signal,length(signal)/nTx,nTx)';
            %
            % ant_1 = sig(1,:);    % symbol txd from antenna 1
            % ant_2 = sig(2,:);    % symbol txd from antenna 2
            
            ant11_conv = conv(signal_1,h11);   % convolution of each symbol with channel h11
            ant11_conv = ant11_conv(1:length(signal_1));
            
            ant12_conv = conv(signal_2,h12);   % convolution of each symbol with channel h12
            ant12_conv = ant12_conv(1:length(signal_2));
            
            ant21_conv = conv(signal_1,h21);   % convolution of each symbol with channel h21
            ant21_conv = ant21_conv(1:length(signal_1));
            
            ant22_conv = conv(signal_2,h22);   % convolution of each symbol with channel h22
            ant22_conv = ant22_conv(1:length(signal_2));
            
            n1 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 1
            n2 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 2
            
            y1 = ant11_conv + ant12_conv + 10^(-SNR(k)/20)*n1 ;
            y2 = ant21_conv + ant22_conv + 10^(-SNR(k)/20)*n2 ;
            
            %  y1 = ant11_conv ;
            %  y2 = ant21_conv + ant22_conv ;
            
            y1 = y1(1:deci:end);
            y2 = y2(1:deci:end);
            
            t = (0:length(y1)-1)/Fs_cor;
            
            tx_sig_1 = y1;  % at antenna 1
            
            tx_sig_2 = y2;  % at antenna 2
            
            
        case '8PSK'
            xa=randi(1,L,8);  % 8 PSK
            
            sig = reshape(xa,[],nTx)';
            ant_1 = sig(1,:);    % symbol txd from antenna 1
            
            const=cos((0:7)*2*pi/8)+sqrt(-1)*sin((0:7)*2*pi/8); % 8psk constellation
            sig_const1=const(ant_1+1);
            
            x_real=real(sig_const1);
            x_img=imag(sig_const1);
            
            x_filtered_r = rcosflt(x_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r1 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';
            
            x_filtered_i = rcosflt(x_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i1 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';
            
            signal_I_Q_1 = x_filtered_r1 + 1j*x_filtered_i1;
            
            ant_2 = sig(2,:);    % symbol txd from antenna 2
            
            const=cos((0:7)*2*pi/8)+sqrt(-1)*sin((0:7)*2*pi/8); % 8psk constellation
            sig_const2=const(ant_2+1);
            
            x_real=real(sig_const2);
            x_img=imag(sig_const2);
            
            x_filtered_r = rcosflt(x_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r2 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';
            
            x_filtered_i = rcosflt(x_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i2 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';
            
            signal_I_Q_2 = x_filtered_r2 + 1j*x_filtered_i2;
            
            ant11_conv = conv(signal_I_Q_1,h11);   % convolution of each symbol with channel h11
            ant11_conv = ant11_conv(1:length(signal_I_Q_1));
            
            ant12_conv = conv(signal_I_Q_2,h12);   % convolution of each symbol with channel h12
            ant12_conv = ant12_conv(1:length(signal_I_Q_2));
            
            ant21_conv = conv(signal_I_Q_1,h21);   % convolution of each symbol with channel h21
            ant21_conv = ant21_conv(1:length(signal_I_Q_1));
            
            ant22_conv = conv(signal_I_Q_2,h22);   % convolution of each symbol with channel h22
            ant22_conv = ant22_conv(1:length(signal_I_Q_2));
            
            n1 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 1
            n2 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 2
            
            y1 = ant11_conv + ant12_conv + 10^(-SNR(k)/20)*n1 ;
            y2 = ant21_conv + ant22_conv + 10^(-SNR(k)/20)*n2 ;
            
            %  y1 = ant11_conv ;
            %  y2 = ant21_conv + ant22_conv ;
            
            y1 = y1(1:deci:end);
            y2 = y2(1:deci:end);
            
            tx_sig_1 = y1;  % at antenna 1
            
            tx_sig_2 = y2;  % at antenna 2
            
            
        case 'PI/4QPSK'
            
            xa = randi(1,L,4);          % PI/4 QPSK
            
            sig = reshape(xa,[],nTx)';
            
            ant_1 = sig(1,:);    % symbol txd from antenna 1
            
            const1=cos(0:pi/2:3*pi/2)+sqrt(-1)*sin(0:pi/2:3*pi/2);    % first constellation
            const2=cos(pi/4:pi/2:2*pi)+sqrt(-1)*sin(pi/4:pi/2:2*pi);  % seecond constellation
            
            sig_const=zeros(1,length(ant_1));
            sig_const(1:2:end)=const1(ant_1(1:2:end)+1);
            sig_const(2:2:end)=const2(ant_1(2:2:end)+1);
            
            x_re=real(sig_const);
            x_img=imag(sig_const);
            
            x_filtered_r_pi = rcosflt(x_re,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r1_pi = x_filtered_r_pi(oversamp*delay+1:end-oversamp*delay)';
            
            x_filtered_i_pi = rcosflt(x_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i1_pi = x_filtered_i_pi(oversamp*delay+1:end-oversamp*delay)';
            
            signal_I_Q_1 = x_filtered_r1_pi + 1j*x_filtered_i1_pi;
            
            ant_2 = sig(2,:);    % symbol txd from antenna 2
            
            const1=cos(0:pi/2:3*pi/2)+sqrt(-1)*sin(0:pi/2:3*pi/2);    % first constellation
            const2=cos(pi/4:pi/2:2*pi)+sqrt(-1)*sin(pi/4:pi/2:2*pi);  % seecond constellation
            
            sig_const=zeros(1,length(ant_2));
            sig_const(1:2:end)=const1(ant_2(1:2:end)+1);
            sig_const(2:2:end)=const2(ant_2(2:2:end)+1);
            
            x_re=real(sig_const);
            x_img=imag(sig_const);
            
            x_filtered_r_pi = rcosflt(x_re,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r1_pi = x_filtered_r_pi(oversamp*delay+1:end-oversamp*delay)';
            
            x_filtered_i_pi = rcosflt(x_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i1_pi = x_filtered_i_pi(oversamp*delay+1:end-oversamp*delay)';
            
            signal_I_Q_2 = x_filtered_r1_pi + 1j*x_filtered_i1_pi;
            
            
            % sig = reshape(signal_I_Q,length(signal_I_Q)/nTx,nTx)';
            %
            % ant_1 = sig(1,:);    % symbol txd from antenna 1
            % ant_2 = sig(2,:);    % symbol txd from antenna 2
            
            ant11_conv = conv(signal_I_Q_1,h11);   % convolution of each symbol with channel h11
            ant11_conv = ant11_conv(1:length(signal_I_Q_1));
            
            ant12_conv = conv(signal_I_Q_2,h12);   % convolution of each symbol with channel h12
            ant12_conv = ant12_conv(1:length(signal_I_Q_2));
            
            ant21_conv = conv(signal_I_Q_1,h21);   % convolution of each symbol with channel h21
            ant21_conv = ant21_conv(1:length(signal_I_Q_1));
            
            ant22_conv = conv(signal_I_Q_2,h22);   % convolution of each symbol with channel h22
            ant22_conv = ant22_conv(1:length(signal_I_Q_2));
            
            n1 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 1
            n2 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 2
            
            y1 = ant11_conv + ant12_conv + 10^(-SNR(k)/20)*n1 ;
            y2 = ant21_conv + ant22_conv + 10^(-SNR(k)/20)*n2 ;
            
            %  y1 = ant11_conv ;
            %  y2 = ant21_conv + ant22_conv ;
            
            y1 = y1(1:deci:end);
            y2 = y2(1:deci:end);
            
            tx_sig_1 = y1;  % at antenna 1
            
            tx_sig_2 = y2;  % at antenna 2
            
            
        case 'MSK'
            
            
            xa=2*randi(1,2*L)-1;  % BPSK signal generation
            
            sig = reshape(xa,[],nTx)';
            ant_1 = sig(1,:);    % symbol txd from antenna 1
            
            e_xa=ant_1(1:2:end);
            o_xa=ant_1(2:2:end);
            
            e_filtered_m = rcosflt(e_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            e_filtered_m = e_filtered_m(oversamp*delay+1:end-oversamp*delay)';
            x_shift_r =[e_filtered_m, zeros(1, floor(oversamp/2))];
            
            o_filtered_m = rcosflt(o_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            o_filtered_m = o_filtered_m(oversamp*delay+1:end-oversamp*delay)';
            x_shift_i = [zeros(1, floor(oversamp/2)), o_filtered_m];
            
            signal_1 = x_shift_r + 1j*x_shift_i;
            
            ant_2 = sig(2,:);    % symbol txd from antenna 1
            
            e_xa=ant_2(1:2:end);
            o_xa=ant_2(2:2:end);
            
            e_filtered_m = rcosflt(e_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            e_filtered_m = e_filtered_m(oversamp*delay+1:end-oversamp*delay)';
            x_shift_r =[e_filtered_m, zeros(1, floor(oversamp/2))];
            
            o_filtered_m = rcosflt(o_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            o_filtered_m = o_filtered_m(oversamp*delay+1:end-oversamp*delay)';
            x_shift_i = [zeros(1, floor(oversamp/2)), o_filtered_m];
            
            signal_2 = x_shift_r + 1j*x_shift_i;
            
            ant11_conv = conv(signal_1,h11);   % convolution of each symbol with channel h11
            ant11_conv = ant11_conv(1:length(signal_1));
            
            ant12_conv = conv(signal_2,h12);   % convolution of each symbol with channel h12
            ant12_conv = ant12_conv(1:length(signal_2));
            
            ant21_conv = conv(signal_1,h21);   % convolution of each symbol with channel h21
            ant21_conv = ant21_conv(1:length(signal_1));
            
            ant22_conv = conv(signal_2,h22);   % convolution of each symbol with channel h22
            ant22_conv = ant22_conv(1:length(signal_2));
            
            n1 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 1
            n2 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 2
            
            y1 = ant11_conv + ant12_conv + 10^(-SNR(k)/20)*n1 ;
            y2 = ant21_conv + ant22_conv + 10^(-SNR(k)/20)*n2 ;
            
            y1 = y1(1:deci:end);
            y2 = y2(1:deci:end);
            
            t = (0:length(y1)-1)/Fs_cor;
            arg = pi*(0:length(y1)-1)/(oversamp/deci);
            
            tx_sig_1 = real(y1).*cos(arg) + 1i*imag(y1).*sin(arg);  % at antenna 1
            
            tx_sig_2 = real(y2).*cos(arg) + 1i*imag(y2).*sin(arg);  % at antenna 2
            
        case '16QAM'
            x_real = 2*randi(1,L,4)-1;
            sig_real = reshape(x_real,[],nTx)';
            ant_1_real = sig_real(1,:);    % symbol txd from antenna 1
            ant_1_real(find(ant_1_real==5)) = -3;
            
            x_img = 2*randi(1,L,4)-1;
            sig_img = reshape(x_img,[],nTx)';
            ant_1_img = sig_img(1,:);    % symbol txd from antenna 1
            ant_1_img(find(ant_1_img==5)) = -3;
            
            
            x_filtered_r = rcosflt(ant_1_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r1 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            x_filtered_i = rcosflt(ant_1_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i1 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            signal_1 = x_filtered_r1 + 1j*x_filtered_i1;
            
            ant_2_real = sig_real(2,:);    % symbol txd from antenna 2
            ant_2_real(find(ant_2_real==5)) = -3;
            
            ant_2_img = sig_img(2,:);    % symbol txd from antenna 1
            ant_2_img(find(ant_2_img==5)) = -3;
            
            x_filtered_r = rcosflt(ant_2_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_r1 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            x_filtered_i = rcosflt(ant_2_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
            x_filtered_i1 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
            
            signal_2 = x_filtered_r1 + 1j*x_filtered_i1;
            
            
            ant11_conv = conv(signal_1,h11);   % convolution of each symbol with channel h11
            ant11_conv = ant11_conv(1:length(signal_1));
            
            ant12_conv = conv(signal_2,h12);   % convolution of each symbol with channel h12
            ant12_conv = ant12_conv(1:length(signal_2));
            
            ant21_conv = conv(signal_1,h21);   % convolution of each symbol with channel h21
            ant21_conv = ant21_conv(1:length(signal_1));
            
            ant22_conv = conv(signal_2,h22);   % convolution of each symbol with channel h22
            ant22_conv = ant22_conv(1:length(signal_2));
            
            n1 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 1
            n2 = 1/sqrt(2)*[randn(size(ant22_conv)) + 1j*randn(size(ant22_conv))]; % white gaussian noise, 0dB variance at 2
            
            y1 = ant11_conv + ant12_conv + 10^(-SNR(k)/20)*n1 ;
            y2 = ant21_conv + ant22_conv + 10^(-SNR(k)/20)*n2 ;
            
            norm_factor = 3/(2*(16-1));
            
            y1 = y1(1:deci:end);
            y2 = y2(1:deci:end);
            
            tx_sig_1 = sqrt(norm_factor/2)*y1;  % at antenna 1
            
            tx_sig_2 = sqrt(norm_factor/2)*y2;  % at antenna 2
            
        otherwise
    end
end
end

