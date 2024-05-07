function [tx_sig] = generate_transmit_signal_fading(L,fc,oversamp,deci,Fs_cor,Mary,fe,theta,R,SNR)


for k=1:length(SNR)
%fadign channel

        sym_block=L;
        n1 = randn(1, L/sym_block)/sqrt(2);% number of symbol where fading will be the same 
        n2 = randn(1, L/sym_block)/sqrt(2);        
        n=sqrt(n1.^2+n2.^2);
        r=kron(n, ones(1, sym_block*oversamp));
        delay = 3;


rc = rcosine(1,oversamp,'sqrt',R);
    
switch Mary

    case 'BPSK'
      
        x_real = 2*randi(1,L)-1;
        x_filtered = rcosflt(x_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        x_filtered1 = x_filtered(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
        x_filtered1_fd= r.*x_filtered1;
        
        x_deci= x_filtered1_fd(1:deci:end);   
        
        noise = 1/sqrt(2)* (10^(-SNR(k)/20))*randn(size(x_deci));
        
        t = (0:length(x_deci)-1)/Fs_cor;
        
         
        tx_sig = (x_deci + noise).* cos(2*pi*fc*t + theta);

    case 'QPSK'
         xa=2*randi(1,2*L)-1;  
        x_evn=xa(1:2:end); 
        x_odd=xa(2:2:end); 

        x_filtered_r = rcosflt(x_evn,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        x_filtered_r1 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
        x_filtered_r_fd = r.*x_filtered_r1 ;
        
        x_filtered_i = rcosflt(x_odd,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        x_filtered_i1 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
        x_filtered_i_fd = r.*x_filtered_i1 ;
        
         x_deci_r = x_filtered_r_fd(1:deci:end);
         x_deci_i = x_filtered_i_fd(1:deci:end);
         
       noise1 = sqrt(1/2)*(10^(-SNR(k)/20))*randn(size( x_deci_r));
       noise2 = sqrt(1/2)*(10^(-SNR(k)/20))*randn(size( x_deci_r));
        t = (0:length( x_deci_r)-1)/Fs_cor;

       tx_sig = ( x_deci_r + noise1).* cos(2*pi*fc*t+ theta) +( x_deci_i + noise2).* sin(2*pi*fc*t+ theta);
%         tx_sig = sqrt(1/2)*(x_filtered_r_fd ).* cos(2*pi*fc*t+ theta) + sqrt(1/2)*(x_filtered_i_fd ).* sin(2*pi*fc*t+ theta);
        
    case 'OQPSK'
       
        xa=2*randi(1,2*L)-1;  
        e_xa=xa(1:2:end); 
        o_xa=xa(2:2:end); 

      
        
        e_filtered = rcosflt(e_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        e_filtered1 = e_filtered(oversamp*delay+1:end-oversamp*delay)';
        e_filtered_fd= r.*e_filtered1;
        
        o_filtered = rcosflt(o_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        o_filtered1 = o_filtered(oversamp*delay+1:end-oversamp*delay)';
        o_filtered_fd= r.*o_filtered1;
        
        x_shift_r=[e_filtered_fd, zeros(1, floor(oversamp/2))];
        x_shift_i=[zeros(1, floor(oversamp/2)), o_filtered_fd];
        
        x_deci_e = x_shift_r(1:deci:end);
        x_deci_o = x_shift_i(1:deci:end);
        
        noise11 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size( x_deci_e));
        noise21 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size( x_deci_o));
        
        t = (0:length(x_deci_e)-1)/Fs_cor;
        x_deci_rr =(x_deci_e + noise11).* cos(2*pi*fc*t + theta) ;
        x_deci_ii=(x_deci_o + noise21).* sin(2*pi*fc*t + theta);
            
        tx_sig=x_deci_rr+x_deci_ii;
        

    case 'PI/4QPSK'
        xa=randi(1,L,4);  
        const1=cos(0:pi/2:3*pi/2)+sqrt(-1)*sin(0:pi/2:3*pi/2);    % first constellation
        const2=cos(pi/4:pi/2:2*pi)+sqrt(-1)*sin(pi/4:pi/2:2*pi);  % seecond constellation           

        sig_const=zeros(1,length(xa));
        sig_const(1:2:end)=const1(xa(1:2:end)+1);
        sig_const(2:2:end)=const2(xa(2:2:end)+1);

        x_real=real(sig_const);
        x_img=imag(sig_const);
        
        x_filtered_r_pi = rcosflt(x_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        x_filtered_r1_pi = x_filtered_r_pi(oversamp*delay+1:end-oversamp*delay)';
        x_filtered_rpi_fd = r.*x_filtered_r1_pi ;
        
        x_filtered_i_pi = rcosflt(x_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        x_filtered_i1_pi = x_filtered_i_pi(oversamp*delay+1:end-oversamp*delay)';
        x_filtered_ipi_fd = r.*x_filtered_i1_pi ;
        
         x_deci_rp =  x_filtered_rpi_fd(1:deci:end);
         x_deci_ip = x_filtered_ipi_fd(1:deci:end);
        

        noise12 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size( x_deci_rp));
        noise22 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size(x_deci_ip));
        

        t = (0:length(x_deci_ip)-1)/Fs_cor;

        tx_sig = (x_deci_rp + noise12).* cos(2*pi*fc*t + theta) + ( x_deci_ip + noise22).* sin(2*pi*fc*t + theta);
        

%     case '8PSK'
%         xa=randint(1,L,8);  % 8 PSK
% 
%         const=cos((0:7)*2*pi/8)+sqrt(-1)*sin((0:7)*2*pi/8); % 8psk constellation
%         sig_const=const(xa+1);
% 
%         x_real=real(sig_const);
%         x_img=imag(sig_const);
% 
%         x_filtered_r = rcosflt(x_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
%         x_filtered_r1 = x_filtered_r(oversamp*delay+1:end-oversamp*delay)';
%         x_filtered_r_fd= r.*x_filtered_r1;
%         
%         x_filtered_i = rcosflt(x_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
%         x_filtered_i1 = x_filtered_i(oversamp*delay+1:end-oversamp*delay)';
%         x_filtered_i_fd= r.*x_filtered_i1;
% 
%         noise00 = (10^(-SNR(k)/20))*randn(size( x_filtered_r_fd));
%         
%         noise01 = (10^(-SNR(k)/20))*randn(size( x_filtered_i_fd));
%         
% 
%         t = (0:length( x_filtered_i_fd)-1)/Fs_cor;
% 
%         tx_sig = (x_filtered_r_fd+noise00) .* cos(2*pi*fc*t + theta) + ( x_filtered_i_fd+noise01) .* sin(2*pi*fc*t + theta);

    case 'MSK'
        
      
       xa=2*randi(1,2*L)-1;  % BPSK signal generation 

        e_xa=xa(1:2:end);
        o_xa=xa(2:2:end); 
        
         e_filtered_m = rcosflt(e_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        e_filtered1_m = e_filtered_m(oversamp*delay+1:end-oversamp*delay)';
        e_filtered_fdm= r.*e_filtered1_m;
        
        o_filtered_m = rcosflt(o_xa,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        o_filtered1_m = o_filtered_m(oversamp*delay+1:end-oversamp*delay)';
        o_filtered_fdm= r.* o_filtered1_m;
        
        x_shift_rr=[e_filtered_fdm, zeros(1, floor(oversamp/2))];
        x_shift_ii=[zeros(1, floor(oversamp/2)), o_filtered_fdm];
        
         x_deci_rr =  x_shift_rr(1:deci:end);
         x_deci_ii = x_shift_ii(1:deci:end);
        
        noise13 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size( x_deci_rr));
        noise23 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size(x_deci_rr));
         
      arg = pi*(0:length(x_deci_rr)-1)/(oversamp/deci);
      t = (0:length(x_deci_rr)-1)/Fs_cor;
      
      x_deci_rr1 = (x_deci_rr + noise13).*cos(arg).* cos(2*pi*fc*t) ;
      x_deci_ii1= (x_deci_ii + noise23).*sin(arg).* sin(2*pi*fc*t );
      
      
       
      tx_sig = x_deci_rr1+x_deci_ii1;     
        
    case '16QAM'
        x_real = 2*randi(1,L,4)-1;
        x_real(find(x_real==5)) = -3;
        
        x_img = 2*randi(1,L,4)-1;
        x_img(find(x_img==5)) = -3;
        
        x_filtered_r_16 = rcosflt(x_real,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        x_filtered_r1_16 = x_filtered_r_16(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
        x_filtered_r16_fd = r.*x_filtered_r1_16;
        
        
        x_filtered_i_16 = rcosflt(x_img,1,oversamp,'filter',rc)*sqrt(oversamp);            %filter
        x_filtered_i1_16 = x_filtered_i_16(oversamp*delay+1:end-oversamp*delay)';  %remove filter transients
        x_filtered_i16_fd = r.*x_filtered_i1_16;
        
         x_deci_r16 =  x_filtered_r16_fd(1:deci:end);
         x_deci_i16 = x_filtered_i16_fd(1:deci:end);
         
        noise14 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size( x_deci_i16));
        noise24 = 1/sqrt(2)*(10^(-SNR(k)/20))*randn(size( x_deci_i16));
        
        t = (0:length( x_deci_i16)-1)/Fs_cor;

        norm_factor = 3/(2*(16-1));
        
        tx_sig = sqrt(norm_factor/2)*( x_deci_r16 + noise14).* cos(2*pi*fc*t + theta) + sqrt(norm_factor/2)*( x_deci_i16 + noise24).* sin(2*pi*fc*t + theta);
        
        
    otherwise 
end
end
end

