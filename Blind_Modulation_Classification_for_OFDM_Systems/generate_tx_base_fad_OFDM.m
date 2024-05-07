function [tx_sig_base_ofdm] = generate_tx_base_fad_OFDM(K,n_symbol,Ns,oversamp,Mary,SNR,h11,Nuse,Ncp)


for k=1:length(SNR)

switch Mary

    case 'BPSK'
        x_real = 2*randi([0,1],1,K*n_symbol)-1; 
 x_real_par = reshape(x_real,K,n_symbol).'; 
    
    
b = [x_real_par(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), x_real_par(:,K/2+1:K)];

tim_sig = K*oversamp*ifft(b, Nuse,2);   % taking the IFFT
        
OFDM_tx = [tim_sig(:,end-Ncp+1:end) tim_sig];  % add cyclic prefix 
   
OFDM_tx_ser = reshape(OFDM_tx.', 1,n_symbol*Ns);
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  
        
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = OFDM_conv + 10^(-SNR(k)/20)*n1;

    case 'QPSK'
        
    
xa=(2*randi([0,1],1,K*n_symbol)-1) + 1i*(2*randi([0,1],1,K*n_symbol)-1);   
x_real_par = reshape(xa,K,n_symbol).';
      
b = [x_real_par(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), x_real_par(:,K/2+1:K)];

tim_sig = K*oversamp*ifft(b, Nuse,2);   % taking the IFFT
        
OFDM_tx = [tim_sig(:,end-Ncp+1:end) tim_sig];  % add cyclic prefix 
   
OFDM_tx_ser = reshape(OFDM_tx.', 1,n_symbol*Ns);
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  
        
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = OFDM_conv + 10^(-SNR(k)/20)*n1;
        
    case 'OQPSK'
    
  xa=2*randi([0,1],1,2*K*n_symbol)-1;  
  x_real_par = reshape(xa,2*K,n_symbol).';

        e_xa = x_real_par(:,1:2:2*K); 
        o_xa = x_real_par(:,2:2:2*K); 
        
b_e = [e_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), e_xa(:,K/2+1:K)];  

tim_sige = K*oversamp*ifft(b_e, Nuse,2);   % taking the IFFT
x_shift_r =[tim_sige, zeros(n_symbol, floor(oversamp/2))];
           

b_o = [o_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), o_xa(:,K/2+1:K)];  

tim_sigo = K*oversamp*ifft(b_o, Nuse,2);   % taking the IFFTx_shift_i =[zeros(1, floor(oversamp/2)), tim_sigo];
x_shift_i =[zeros(n_symbol, floor(oversamp/2)) tim_sigo];
      
signal_1 = x_shift_r + 1j*x_shift_i;
             
OFDM_tx = [signal_1(:,end-Ncp+1:end) signal_1];  % add cyclic prefix 

OFDM_tx_ser = reshape(OFDM_tx.', 1,size(OFDM_tx,1)* size(OFDM_tx,2));
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  
        
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = OFDM_conv + 10^(-SNR(k)/20)*n1;

    case 'PI/4QPSK'
          
            
       xa = randi([0,3],1,K*n_symbol);       % PI/4 QPSK 
        x_real_par = reshape(xa,K,n_symbol).';        

        const1=cos(0:pi/2:3*pi/2)+sqrt(-1)*sin(0:pi/2:3*pi/2);    % first constellation
        const2=cos(pi/4:pi/2:2*pi)+sqrt(-1)*sin(pi/4:pi/2:2*pi);  % seecond constellation           

        sig_const=zeros(n_symbol,K);
        sig_const(:,1:2:end)=const1(x_real_par(:,1:2:end)+1);
        sig_const(:,2:2:end)=const2(x_real_par(:,2:2:end)+1);

        e_xa=real(sig_const);
        o_xa=imag(sig_const);
   
b_e = [e_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), e_xa(:,K/2+1:K)];  
tim_sige = K*oversamp*ifft(b_e, Nuse,2);   % taking the IFFT

b_o = [o_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), o_xa(:,K/2+1:K)];  
tim_sigo = K*oversamp*ifft(b_o, Nuse,2);   % taking the IFFTx_shift_i =[zeros(1, floor(oversamp/2)), tim_sigo];

signal_1 = tim_sige + 1j*tim_sigo;
            
OFDM_tx = [signal_1(:,end-Ncp+1:end) signal_1];  % add cyclic prefix 

OFDM_tx_ser = reshape(OFDM_tx.', 1,size(OFDM_tx,1)* size(OFDM_tx,2));
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  
        
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = OFDM_conv + 10^(-SNR(k)/20)*n1;
        
        
        case '8PSK'

        xa=randi([0,7],1,K*n_symbol);  % 8 PSK
        x_real_par = reshape(xa,K,n_symbol).'; 
        
        const=cos((0:7)*2*pi/8)+sqrt(-1)*sin((0:7)*2*pi/8); % 8psk constellation
        sig_const=const(x_real_par+1);

        e_xa=real(sig_const);
        o_xa=imag(sig_const);
   
b_e = [e_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), e_xa(:,K/2+1:K)];  
tim_sige = K*oversamp*ifft(b_e, Nuse,2);   % taking the IFFT

b_o = [o_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), o_xa(:,K/2+1:K)];  
tim_sigo = K*oversamp*ifft(b_o, Nuse,2);   % taking the IFFTx_shift_i =[zeros(1, floor(oversamp/2)), tim_sigo];

signal_1 = tim_sige + 1j*tim_sigo;
            
OFDM_tx = [signal_1(:,end-Ncp+1:end) signal_1];  % add cyclic prefix 

OFDM_tx_ser = reshape(OFDM_tx.', 1,size(OFDM_tx,1)* size(OFDM_tx,2));
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  
        
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = OFDM_conv + 10^(-SNR(k)/20)*n1;

    case 'MSK'
        
  xa=2*randi([0,1],1,2*K*n_symbol)-1;  
  x_real_par = reshape(xa,2*K,n_symbol).';

        e_xa = x_real_par(:,1:2:2*K); 
        o_xa = x_real_par(:,2:2:2*K); 
        
 
b_e = [e_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), e_xa(:,K/2+1:K)];  
tim_sige = K*oversamp*ifft(b_e, Nuse,2);   % taking the IFFT
x_shift_r =[tim_sige, zeros(n_symbol, floor(oversamp/2))];
           

b_o = [o_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), o_xa(:,K/2+1:K)];  
tim_sigo = K*oversamp*ifft(b_o, Nuse,2);   % taking the IFFTx_shift_i =[zeros(1, floor(oversamp/2)), tim_sigo];
x_shift_i =[zeros(n_symbol, floor(oversamp/2)) tim_sigo];
      
signal_1 = x_shift_r + 1j*x_shift_i;
             
OFDM_tx = [signal_1(:,end-Ncp+1:end) signal_1];  % add cyclic prefix 

OFDM_tx_ser = reshape(OFDM_tx.', 1,size(OFDM_tx,1)* size(OFDM_tx,2));
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));            
        
        arg = pi*(0:length(OFDM_conv)-1)/(oversamp);
    
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        x_deci_rr1 = OFDM_conv + 10^(-SNR(k)/20)*n1;     
             
      tx_sig_base_ofdm = real(x_deci_rr1).*cos(arg) + 1i*imag(x_deci_rr1).*sin(arg);  % at antenna 1   
  
    case '16QAM'
          
 
        x_real  = 2*randi([0,3],1,K*n_symbol) -1;        
        x_real(find(x_real==5)) = -3;
        x_real_par = reshape(x_real,K,n_symbol).';

        
        x_img = 2*randi([0,3],1,K*n_symbol) -1;      
        x_img(find(x_img==5)) = -3;
        x_img_par = reshape(x_img,K,n_symbol).';  
          
b_e = [x_real_par(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), x_real_par(:,K/2+1:K)];  
tim_sige = K*oversamp*ifft(b_e, Nuse,2);   % taking the IFFT

b_o = [x_img_par(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), x_img_par(:,K/2+1:K)];  
tim_sigo = K*oversamp*ifft(b_o, Nuse,2);   % taking the IFFTx_shift_i =[zeros(1, floor(oversamp/2)), tim_sigo];

signal_1 = tim_sige + 1j*tim_sigo;
            
OFDM_tx = [signal_1(:,end-Ncp+1:end) signal_1];  % add cyclic prefix 

OFDM_tx_ser = reshape(OFDM_tx.', 1,size(OFDM_tx,1)* size(OFDM_tx,2));
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  

         norm_factor = 3/(2*(16-1));
   
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = sqrt(norm_factor/2)*(OFDM_conv + 10^(-SNR(k)/20)*n1);        
        
         case '64QAM'
       
        data = randi([0,63],n_symbol,K);
        sig_const = (qammod(data,64));    % for normalization divide by sqrt(42)

        e_xa=real(sig_const);
        o_xa=imag(sig_const);
   
b_e = [e_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), e_xa(:,K/2+1:K)];  
tim_sige = K*oversamp*ifft(b_e, Nuse,2);   % taking the IFFT

b_o = [o_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), o_xa(:,K/2+1:K)];  
tim_sigo = K*oversamp*ifft(b_o, Nuse,2);   % taking the IFFTx_shift_i =[zeros(1, floor(oversamp/2)), tim_sigo];

signal_1 = tim_sige + 1j*tim_sigo;
            
OFDM_tx = [signal_1(:,end-Ncp+1:end) signal_1];  % add cyclic prefix 

OFDM_tx_ser = reshape(OFDM_tx.', 1,size(OFDM_tx,1)* size(OFDM_tx,2));
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  
        
        norm_factor = 3/(2*(64-1));

        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = sqrt(norm_factor/2)*(OFDM_conv + 10^(-SNR(k)/20)*n1);   
        
        
         case '16PSK'
             
      xa=randi([0,15],1,K*n_symbol);  % 8 PSK
        x_real_par = reshape(xa,K,n_symbol).'; 
        
        const=cos((0:15)*2*pi/16)+sqrt(-1)*sin((0:15)*2*pi/16); % 8psk constellation
        sig_const=const(x_real_par+1);

        e_xa=real(sig_const);
        o_xa=imag(sig_const);
   
b_e = [e_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), e_xa(:,K/2+1:K)];  
tim_sige = K*oversamp*ifft(b_e, Nuse,2);   % taking the IFFT

b_o = [o_xa(:,1:K/2), zeros(n_symbol, K*(oversamp-1)), o_xa(:,K/2+1:K)];  
tim_sigo = K*oversamp*ifft(b_o, Nuse,2);   % taking the IFFTx_shift_i =[zeros(1, floor(oversamp/2)), tim_sigo];

signal_1 = tim_sige + 1j*tim_sigo;
            
OFDM_tx = [signal_1(:,end-Ncp+1:end) signal_1];  % add cyclic prefix 

OFDM_tx_ser = reshape(OFDM_tx.', 1,size(OFDM_tx,1)* size(OFDM_tx,2));
       
OFDM_conv = conv(OFDM_tx_ser,h11);            % convolution of each symbol with channel h
OFDM_conv = OFDM_conv(1:length(OFDM_tx_ser));  
        
        n1 = 1/sqrt(2)*[randn(size(OFDM_conv)) + 1j*randn(size(OFDM_conv))]; % white gaussian noise, 0dB variance at 1
        tx_sig_base_ofdm = OFDM_conv + 10^(-SNR(k)/20)*n1;   
        
    otherwise 
end
end
end

