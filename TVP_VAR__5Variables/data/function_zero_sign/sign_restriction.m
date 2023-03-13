function [P,FP,MP,BC,chk] = sign_restriction(amOmsq,nk, t, nl,nlen,mbs)

chk =1;

% nk = 系列の数
% amOmsq = inv(A)*Sigma = 1期目のレスポンス応答
% t = ショックをおこすperiod　(nl+1(ラグの次数+1)～ns(データのサンプル期間)) 

a=0; % Nov30, 2013
b=0; % Nov30, 2013
d=0;

my = zeros(nl+nlen, nk);   % インパルスレスポンス

while d==0

% Omega = QR decompotion
% RWZ (2010) p688, algorithm 2 

x = randn(nk,nk);

[P, R]=qr(x);        %  i.e., P * P' = I

i = 1;    % 財政政策(政府ショック)ショックを起こす係数のorder 
%  my(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1期目の制約付きレスポンス
my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;

%  2期以降の　インパルス応答
  for j = nl+2 : nl+nlen
       my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
  end
  
  Imp_FP = zeros(nk,nlen);
  Imp_FP =  my(nl+1: nl+nlen,:)' ;

%  系列 1 = gov cons & inv
%  系列 2 = real GDP
%  系列 3 = private consumption
%  系列 4 = Debt GDP ratio
%  系列 5 = GDP deflator
%  系列 6 = Interest Rate
%  系列 7 = real exchange rate
%  系列 8 = net export

%  sign_restrictionの設定
%  Imp_t1(系列数)

if( Imp_FP(1,1:1)>0) &(Imp_FP(2,1:1)>0) &(Imp_FP(4,1:1)>0);
    FP = 1;     a=1; % Nov30, 2013
elseif (-1*Imp_FP(1,1:1)>0 )&(-1*Imp_FP(2,1:1)>0)&(-1*Imp_FP(4,1:1)>0); % Dec 11, 2013
    FP = -1;    a=1; % Dec 11, 2013

    
    % 財政政策ショックに対して政府支出(系列数;1)の反応は 1期プラス
    % 財政政策ショックに対してGDP(系列数;2)の反応は 1期プラス
    % 財政政策ショックに対してDEBT(系列数;4)の反応は 1期プラス   
    
else
    a=0; % Nov30, 2013
end

% Nov30, 2013 %%%%%%%%%%

% a=1; FP = 1;

if a==0
    d=0;
    
else     
    
i = 3;    % 金融政策(利子率ショック)ショックを起こす係数のorder 
%  my2(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1期目の制約付きレスポンス
 my2(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;

%  2期以降の　インパルス応答
  for j = nl+2 : nl+nlen
       my2(j, :) = mbs(t+j-nl-1,:) * fXt(my2(j-nl:j-1,:), 0)';
  end
  
  Imp_MP = zeros(nk,nlen);
  Imp_MP =  my2(nl+1: nl+nlen,:)' ;

    
 if (Imp_MP(2,1:1)<0 )&( Imp_MP(5,1:1)<0 )&( Imp_MP(6,1:1)>0 );
                     MP = 1;      b=1; 
elseif (-1*Imp_MP(2,1:1)<0)& (-1*Imp_MP(5,1:1)<0)& (-1*Imp_MP(6,1:1)>0);
      MP = -1;     b=1;    

    % 金融引締めショックに対してGDP(系列数;2)の反応は 1期マイナス
    % 金融引締めショックに対してインフレ率(系列数;5)の反応は 1期マイナス
    % 金融引締めショックに対して利子率(系列数;6)の反応は 1期プラス
  
else
    b=0;
 end
 
% Dec13, 2013 %%%%%%%%%%

% a=1; FP = 1;

if b==0
    d=0;
    
else     
    
i = 2;    % 生産性ショックを起こす係数のorder 
%  my3(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1期目の制約付きレスポンス
 my3(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;

%  2期以降の　インパルス応答
  for j = nl+2 : nl+nlen
       my3(j, :) = mbs(t+j-nl-1,:) * fXt(my3(j-nl:j-1,:), 0)';
  end
  
  Imp_BC = zeros(nk,nlen);
  Imp_BC =  my3(nl+1: nl+nlen,:)' ;

if (Imp_BC(2,1:1)>0 )&( Imp_BC(3,1:1)>0 )&( Imp_BC(4,1:1)<0 );
                     BC = 1;      d=1; 
elseif (-1*Imp_BC(2,1:1)>0)&(-1*Imp_BC(3,1:1)>0 )& (-1*Imp_BC(4,1:1)<0);
      BC = -1;     d=1;    
    % 生産性ショックに対してGDP(系列数;2)の反応は 1期プラス
    % 生産性ショックに対して消費(系列数;3)の反応は 1期プラス
    % 生産性ショックに対してDEBT(系列数;4)の反応は 1期マイナス
  
else
    d=0;
 end
end
%%%%%%%%%%%%%%%%%%%%%%
      
end
end
% P_FP = FP*P;
% 
% P_MP = MP*P;

