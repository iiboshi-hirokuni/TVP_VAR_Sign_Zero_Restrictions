function [P, sign_P, sign_chk,A0_t] = ...
    zero_sign_restriction_v4(amOmsq,nk, t, nl,nlen,mbs,svar_type)

  sign_P=ones(nk,1);
  
%====================================================
%  Make Matrix Q and S for Zero and Sign Restrictions
%====================================================
           P=[];
           IR_0=[];
           IR_Long=[];
           Phi0= reshape(mbs(t,:)',nk*nl,nk); P0 = eye(nk);
           
       for i = 1:nk 
          my(nl+1, :) = amOmsq(:,i,t)' ;  % 1期目の制約付きレスポンス          
         
        %  2期以降の　インパルス応答
           for j = nl+2 : nl+nlen       
%                h = j-nl; 
%                my(j,:) = impulseresponce(Phi0,amOmsq(:,:,t),P0(:,i),h);  
                  my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
           end          
  
         IR_0 = [IR_0;  my(nl+1, :) ];       % 短期のインパルス応答
         IR_Long = [IR_Long; my(nl+nlen,:)];  % 長期のインパルス応答 (RWT,2010, P670)
         
       end
      %disp('短期と長期のインパルス応答　関数 F');
%            IR_0= amOmsq(:,:,t)';         
%          n = size(Phi0,2);
%          A_L =zeros(n,n);
%          A0 = inv(amOmsq(:,:,t)); 
%          
%          for i = 1:nl
%            B= Phi0(1+n*(i-1):n*i,:);   
%            A_L = A_L + B*A0;
%          end   
% 
%         % 長期のインパルス応答 (RWT,2010, P670)
%          IR_inf = inv(A0' - A_L');
 
    A0 = inv(amOmsq(:,:,t)'); 
      
switch svar_type
    case 'Caldala_Kamp'
      definition_matrix_Q_CK;      
     % 行列: A_0 の制約
%        A0 = inv(amOmsq(:,:,t)'); 
      F = [A0; IR_0'; IR_Long' ];
%        F = [A0; IR_Long' ];

   case 'Caldala_Kamps_tax'
      definition_matrix_Q_CK_tax;      
     % 行列: A_0 の制約
%        A0 = inv(amOmsq(:,:,t)'); 
      F = [A0; IR_0'; IR_Long' ];
%        F = [A0; IR_Long' ];
 
   case 'Mountford_Uhlig'  
        definition_matrix_Q_MU; 
 %   A0 = inv(amOmsq(:,:,t)');  
       F = [IR_0'; IR_Long'];
%        F = [IR_0'; IR_inf'];
%      xx = randn(nk,nk);    
%     [P, T] = qr(xx);
   end

%====================================================
%  Zero Restriction
%====================================================
% 
% switch svar_type
%  case 'Caldala_Kamp'

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
   
%   case 'Mountford_Uhlig' 
% 
%  for i = 1:(nk-1)
% 
%    Q_bar = abs(orth(squeeze(Q(:,:,i))')); 
%    
%    if isempty(Q_bar) == 1
%          x1 = randn(nk-i,nk); % 制約iをダミーで処理
%          Q_til = [x1; P'];
%          [QQ,T]=qr(Q_til');
%    else    
%          % 制約 i の処理 
%            Q_til = [Q_bar' * F; P'];
%          % QR decomposition
%            [QQ,T]=qr(Q_til');
%    end
%          col_P = QQ(:,nk); % 直交行列 Pのi列目
%          P=[P col_P ]; 
%  end    
%    
%    Q_til = P';
%    [QQ,T]=qr(Q_til');
%    col_P = QQ(:,nk); % 直交行列 Pの最終列目
%    P=[P col_P ];  
% end
       
    
%====================================================
%  Sign Restriction
%====================================================

sign_OK = 1; 
A0_t = zeros(nk,1);

switch svar_type
 case 'Caldala_Kamp'
  Sign_Res_CK;    
  
 case 'Caldala_Kamps_tax'
  Sign_Res_CK_tax;   
  
case 'Mountford_Uhlig'  
  Sign_Res_MU;
end  
    
if sign_OK == 1
     sign_chk = 1; % Sign Restriction Pass
else
     sign_chk = 0;  % Sign Restriction fail    
end

end % end of function
