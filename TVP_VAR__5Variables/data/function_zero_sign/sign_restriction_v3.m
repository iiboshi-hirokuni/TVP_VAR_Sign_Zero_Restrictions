function [P, sign_P, sign_chk] = sign_restriction_v3(amOmsq,nk, t, nl,nlen,mbs)

  sign_P=ones(nk,1);
  
%====================================================
%  Make Matrix Q and S for Zero and Sign Restrictions
%====================================================

definition_matrix_Q_S_v2; 

IR_0=[];
IR_inf=[];

for i = 1:nk 
 my(nl+1, :) = amOmsq(:,i,t)' ;  % 1���ڂ̐���t�����X�|���X

%  2���ȍ~�́@�C���p���X����
  for j = nl+2 : nl+nlen
       my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
  end
  
  IR_0 = [IR_0;  my(nl+1, :) ];       % �Z���̃C���p���X����
  IR_inf = [IR_inf; my(nl+nlen,:)];  % �����̃C���p���X���� (RWT,2010, P670)
end

%disp('�Z���ƒ����̃C���p���X�����@�֐� F');
F = [IR_0'; IR_inf'];

    
%====================================================
%  Sign Restriction
%====================================================

sign_OK = 0; 
%  sign_OK = 1;  FP =1; MP=1; BC=1;
% 

while sign_OK == 0
    
     sign_OK = 0; 
    
      xx = randn(nk,nk);
    
    [P, T] = qr(xx);
    
for i = 1:3 % 
%     S_bar = orth(squeeze(S(:,:,i))');
    
    S = F*P(:,i);
    if i ==1  %  ��������V���b�N
      if (S(1)>0)&&(S(2)>0)&&(S(4)>0)
         sign_P(i) = 1; 
         sign_OK = 1;       
      elseif   (S(1)<0)&&(S(2)<0)&&(S(4)<0)
          sign_P(i) = -1;
          sign_OK = 1;
      else   
       sign_OK = 0; % Sign Restriction fail
        break;
      end
%     end   
    elseif i ==2  %  �i�C�z�V���b�N
      if (S(2)>0)&&(S(3)>0)&&(S(4)<0)
         sign_P(i) = 1; 
         sign_OK = 1;       
      elseif   (S(2)<0)&&(S(3)<0)&&(S(4)>0)
          sign_P(i) = -1;
          sign_OK = 1;
      else   
       sign_OK = 0; % Sign Restriction fail
        break;
      end
%   end  
    else %if i ==3 %  ���Z����V���b�N
      if (S(2)>0)&&(S(5)>0)&&(S(6)<0)
         sign_P(i) = 1; 
         sign_OK = 1;       
      elseif   (S(2)<0)&&(S(5)<0)&&(S(6)>0)
          sign_P(i) = -1;
          sign_OK = 1;
      else   
       sign_OK = 0; % Sign Restriction fail
        break;
      end
   end
   
  end
end


if sign_OK == 1;
     sign_chk = 1; % Sign Restriction Pass
else
     sign_chk = 0;  % Sign Restriction fail    
end

end % end of function
