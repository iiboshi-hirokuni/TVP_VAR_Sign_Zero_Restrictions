
%  Mountford_Uhlig


for i = 1:7 % 2:MP, 3:FP, 4:BC  

%     S_bar = orth(squeeze(S(:,:,i))');
    
    S = F*P(:,i);
    if i == 1   % �i�C�z�V���b�N���N�����W����order 
       if (S(1)>0)&&(S(2)>0)&&(S(3)>0)&&(S(5)>0) % Y >0, C>0, Tax>0, Inv>0     
           % 
         sign_P(i) = 1; 
         sign_OK = 1;       
        elseif   (S(1)<0)&&(S(2)<0)&&(S(3)<0) &&(S(5)<0)  
          sign_P(i) = -1;
          sign_OK = 1;
      else   
       sign_OK = 0; % Sign Restriction fail
        break;
      end
%     end   
    elseif i == 2  % ���Z����V���b�N���N�����W����order 
      if (S(6)>0)&&(S(7)<0) %  % r >0 & Pi < 0
         sign_P(i) = 1; 
         sign_OK = 1;       
      elseif   (S(6)<0)&&(S(7)>0) 
          sign_P(i) = -1;
          sign_OK = 1;
      else   
          sign_OK = 0; % Sign Restriction fail
          break;
      end 
    elseif i == 3 % ��������(���{�V���b�N)�V���b�N���N�����W����order
        if (S(4)>0) % Gov >0
         sign_P(i) = 1; 
         sign_OK = 1;   
           A0_t = A0*P(:,i);
       elseif  (S(4)<0) % 
          sign_P(i) = -1;
          sign_OK = 1;  
          A0_t = -1*A0*P(:,i);
       else   
          sign_OK = 0; % Sign Restriction fail
          break;
        end
    elseif i == 4 % �ŋ��V���b�N���N�����W����order
        if (S(3)>0) % tax >0
         sign_P(i) = 1; 
         sign_OK = 1;       
       elseif  (S(3)<0)  
          sign_P(i) = -1;
          sign_OK = 1;        
       else   
          sign_OK = 0; % Sign Restriction fail
          break;
        end 
        
        
    end 
   
  end % for i = 2:4