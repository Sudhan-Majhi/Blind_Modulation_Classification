function [mod_class] = mod_class_algorithm(oversamp,C_tile_42_d,C_tile_42_fft,C_tile_42_dwt_db) 

if (19 <= oversamp) && (oversamp <= 21)

           if (C_tile_42_d > 33.5)   
               mod_class='OQPSK';                   
           elseif  (C_tile_42_d < 22 ) 
                mod_class='MSK';                            
           elseif   (29.5 <= C_tile_42_d ) && (C_tile_42_d <= 33.5) 
                mod_class='16QAM';                   
           elseif   (22 <= C_tile_42_d ) && (C_tile_42_d < 29.5) &&  (C_tile_42_fft > 34) 
                mod_class='BPSK';             
           elseif   (22 <= C_tile_42_d ) && (C_tile_42_d < 29.5) &&  (C_tile_42_fft <= 34) && (C_tile_42_dwt_db <= 44 )
                mod_class='QPSK';            
           else (22 <= C_tile_42_d ) && (C_tile_42_d < 29.5) &&  (C_tile_42_fft <= 34) && (C_tile_42_dwt_db > 44 );
               mod_class='PI/4QPSK';
           end
           
elseif (24 <= oversamp) && (oversamp <= 26)    
    
           if (C_tile_42_d > 41)   
               mod_class='OQPSK';                   
           elseif  (C_tile_42_d < 28.5 ) 
                mod_class='MSK';                            
           elseif   ( 36.5<= C_tile_42_d ) && (C_tile_42_d <= 41) 
                mod_class='16QAM';                    
           elseif   (28.5 <= C_tile_42_d ) && (C_tile_42_d < 36.5 ) &&  (C_tile_42_fft > 41) 
                mod_class='BPSK';             
           elseif   (28.5 <= C_tile_42_d ) && (C_tile_42_d < 36.5) &&  (C_tile_42_fft <= 41) && (C_tile_42_dwt_db <= 57 )
                mod_class='QPSK';  
           else     (28.5 <= C_tile_42_d ) && (C_tile_42_d < 36.5) &&  (C_tile_42_fft <= 41) && (C_tile_42_dwt_db > 57 );
               mod_class='PI/4QPSK';
           end
           
           
elseif (29 <= oversamp) && (oversamp <= 31)
               
           if (C_tile_42_d > 50)   
               mod_class='OQPSK';                    
           elseif  (C_tile_42_d < 34 ) 
                mod_class='MSK';                            
           elseif   ( 44 < C_tile_42_d ) && (C_tile_42_d <= 50) 
                mod_class='16QAM';                    
           elseif   (34 <= C_tile_42_d ) && (C_tile_42_d <= 44 ) &&  (C_tile_42_fft > 50) 
                mod_class='BPSK';                
           elseif   (34 <= C_tile_42_d ) && (C_tile_42_d <= 44) &&  (C_tile_42_fft <= 50) && (C_tile_42_dwt_db <= 68 )
                mod_class='QPSK';            
           else (34 <= C_tile_42_d ) && (C_tile_42_d <= 44) &&  (C_tile_42_fft <= 50) && (C_tile_42_dwt_db > 68 );
               mod_class='PI/4QPSK';
           end
           
elseif (39 <= oversamp) && (oversamp <= 41)

           if (C_tile_42_d > 67.5)   
               mod_class='OQPSK';                    
           elseif  (C_tile_42_d < 47 ) 
                mod_class='MSK';                           
           elseif   ( 60 < C_tile_42_d ) && (C_tile_42_d <= 67.5) 
                mod_class='16QAM';                     
           elseif   (47 <= C_tile_42_d ) && (C_tile_42_d <= 60 ) &&  (C_tile_42_fft > 66) 
                mod_class='BPSK';                
           elseif   (47 <= C_tile_42_d ) && (C_tile_42_d <= 60) &&  (C_tile_42_fft <= 66) && (C_tile_42_dwt_db <= 90 )
                mod_class='QPSK';           
           else (47 <= C_tile_42_d ) && (C_tile_42_d <= 60) &&  (C_tile_42_fft <= 66) && (C_tile_42_dwt_db > 90 );
               mod_class='PI/4QPSK';
           end         
       
else (49 <= oversamp) && (oversamp <= 51);

           if (C_tile_42_d > 84)   
              mod_class='OQPSK';                    
           elseif  (C_tile_42_d < 60) 
                mod_class='MSK';                           
           elseif   (75 <= C_tile_42_d ) && (C_tile_42_d <= 84) 
                mod_class='16QAM';                     
           elseif   (60 <= C_tile_42_d ) && (C_tile_42_d < 75) &&  (C_tile_42_fft > 82) 
                mod_class='BPSK';               
           elseif   (60 <= C_tile_42_d ) && (C_tile_42_d < 75) &&  (C_tile_42_fft <= 82) && (C_tile_42_dwt_db <= 111 )
                mod_class='QPSK';           
           else (60 <= C_tile_42_d ) && (C_tile_42_d < 75) &&  (C_tile_42_fft <= 82) && (C_tile_42_dwt_db > 111 );
               mod_class='PI/4QPSK';
           end
                   
end
end