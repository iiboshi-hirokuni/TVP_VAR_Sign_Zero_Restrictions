function [P, sign_P, sign_chk] = zero_sign_restriction_v31(amOmsq,nk, t, nl,nlen,mbs)

  sign_P=ones(nk,1);

  global svar_type;
  switch svar_type
      case 'both_zero_sign'    
            type_of_restriction = 'Zero_Restriction';     % �[������ƕ�������̂Q�d������s��
      case 'just_sign'
            type_of_restriction = 'Just_Sign_Restriction';  % ��������̂ݍs��
  end
  
%====================================================
%  Make Matrix Q and S for Zero and Sign Restrictions
%====================================================

definition_matrix_Q_S_v3; 

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

P=[];

%====================================================
%  Zero Restriction
%====================================================

switch type_of_restriction
case 'Zero_Restriction'

 for i = 1:(nk-1)

   Q_bar = abs(orth(squeeze(Q(:,:,i))')); 
   
   if isempty(Q_bar) == 1
         x1 = randn(nk-i,nk); % ����i���_�~�[�ŏ���
         Q_til = [x1; P'];
         [QQ,T]=qr(Q_til');
   else    
         % ���� i �̏��� 
           Q_til = [Q_bar' * F; P'];
         % QR decomposition
           [QQ,T]=qr(Q_til');
   end
         col_P = QQ(:,nk); % �����s�� P��i���
         P=[P col_P ]; 
 end    
   
   Q_til = P';
   [QQ,T]=qr(Q_til');
   col_P = QQ(:,nk); % �����s�� P�̍ŏI���
   P=[P col_P ]; 
   
case 'Just_Sign_Restriction'
     
    xx = randn(nk,nk);
    
    [P, T] = qr(xx);
end
       
    
%====================================================
%  Sign Restriction
%====================================================

sign_OK = 0;  FP =0; MP=0; BC=0;
%  sign_OK = 1;  FP =1; MP=1; BC=1;
% 
for i = 3:5 % 
%     S_bar = orth(squeeze(S(:,:,i))');
    
    S = F*P(:,i);
    if i == 3  %  ��������V���b�N
      if (S(1)>0)&&(S(2)>0)&&(S(4)>0)&&(S(5)>0)&&(S(6)>0)
         sign_P(i) = 1; 
         sign_OK = 1;       
      elseif   (S(1)<0)&&(S(2)<0)&&(S(4)<0)&&(S(5)<0)&&(S(6)<0)
          sign_P(i) = -1;
          sign_OK = 1;
      else   
       sign_OK = 0; % Sign Restriction fail
        break;
      end
%     end   
    elseif i == 4  %  �i�C�z�V���b�N
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
    else %if i == 1 %  ���Z����V���b�N
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

if sign_OK == 1;
     sign_chk = 1; % Sign Restriction Pass
else
     sign_chk = 0;  % Sign Restriction fail    
end

end % end of function
