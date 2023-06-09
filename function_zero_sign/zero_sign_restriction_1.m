function [P, FP, BC, MP, sign_chk] = zero_sign_restriction_1(amOmsq,nk, t, nl,nlen,mbs)

%   sign_P=ones(nk,1);

%  type_of_restriction = 'Zero_Restriction';     % ゼロ制約と符号制約の２重制約を行う

  type_of_restriction = 'Just_Sign_Restriction';  % 符号制約のみ行う

%====================================================
%  Make Matrix Q and S for Zero and Sign Restrictions
%====================================================

definition_matrix_Q_S; 

IR_0=[];
IR_inf=[];

for i = 1:nk 
 my(nl+1, :) = amOmsq(:,i,t)' ;  % 1期目の制約付きレスポンス

%  2期以降の　インパルス応答
  for j = nl+2 : nl+nlen
       my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
  end
  
  IR_0 = [IR_0;  my(nl+1, :) ];       % 短期のインパルス応答
  IR_inf = [IR_inf; my(nl+nlen,:)];  % 長期のインパルス応答 (RWT,2010, P670)
end

%disp('短期と長期のインパルス応答　関数 F');
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
         x1 = randn(nk-i,nk); % 制約iをダミーで処理
         Q_til = [x1; P'];
         [QQ,T]=qr(Q_til');
   else    
         % 制約 i の処理 
           Q_til = [Q_bar' * F; P'];
         % QR decomposition
           [QQ,T]=qr(Q_til');
   end
         col_P = QQ(:,nk); % 直交行列 Pのi列目
         P=[P col_P ]; 
 end    
   
   Q_til = P';
   [QQ,T]=qr(Q_til');
   col_P = QQ(:,nk); % 直交行列 Pの最終列目
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
for i = 1:size(S,3)
    S_bar = orth(squeeze(S(:,:,i))');
    
%   if isempty(S_bar)== 0 % S_barは空ではない   
    Sign_Res = S_bar'*F*P(:,i); 

    if (Sign_Res > 0)
%         sign_P(i) = 1; 
        sign_OK = 1;       
       if i ==1  FP = 1;    
       elseif  i ==2  BC = 1;
       elseif  i ==3  MP = 1;   end
    elseif   (-1*Sign_Res > 0)
%           sign_P(i) = -1;
          sign_OK = 1;
       if i ==1  FP = -1; 
       elseif  i ==2  BC = -1;
       elseif  i ==3  MP = -1;   end
    else   
       sign_OK = 0; % Sign Restriction fail
       break;
    end
%     end
    
end

if sign_OK == 1;
     sign_chk = 1; % Sign Restriction Pass
else
     sign_chk = 0;  % Sign Restriction fail    
end

end % end of function
