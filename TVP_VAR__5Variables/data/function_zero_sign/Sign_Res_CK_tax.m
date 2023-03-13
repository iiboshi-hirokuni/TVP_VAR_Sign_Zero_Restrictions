
%  Caldara and Kamps (2017)

 for i = 1:5 % 2:MP, 3:FP, 4:BC  
    
    S = F*P(:,i);
    if  i == 3 % 財政政策(Taxショック)ショックを起こす係数のorder
        if (S(1)<0)&&(S(2)>0)&&(S(3)>0)&&(S(4)<0)&&(S(5)<0)&&(S(2+5)>0) % y<0 Tax>0, Gov>0 R<0, Pi<0, short Tax>0
            A0_t = S(1:nk);
%            sign_P(i) = 1/A0_t(4);   % omega_g = 1/a0_11 (Caldara and Kamps, 2017, page 1019)
            sign_P(i) = 1;
               
       elseif (S(1)>0)&&(S(2)<0)&&(S(3)<0)&&(S(4)>0)&&(S(5)>0)&&(S(2+5)<0)  % y>0 Tax<0, Gov<0 R>0, Pi>0, short Tax<0
            A0_t = -S(1:nk);  
%           sign_P(i) = -1/A0_t(4);   
           sign_P(i) = -1; 
       else   
          sign_OK = 0; % Sign Restriction fail
          break;
        end

    end 
   
   end % for i = 2:4