function [P] = sign_restriction(amOmsq,nk, t, nl,nlen,mbs)

% nk = 系列の数
% amOmsq = inv(A)*Sigma = 1期目のレスポンス応答
% t = ショックをおこすperiod　(nl+1(ラグの次数+1)～ns(データのサンプル期間)) 

d=0;

my = zeros(nl+nlen, nk);   % インパルスレスポンス

while d==0

%  Omega = QR decompotion
% RWZ (2010) p688, algorithm 2 

x = randn(nk,nk);

[P, R]=qr(x);        %  i.e., P * P' = I

i = 1;    % 財政政策(政府ショック)ショックを起こす係数のorder 
 my(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1期目の制約付きレスポンス

%  2期以降の　インパルス応答
  for j = nl+2 : nl+nlen
       my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
  end
  
  Imp_t = zeros(nk,nlen);
  Imp_t =  my(nl+1: nl+nlen,:)' ;

% 　系列 1 = 政府消費
%   系列 2 = real GDP
%   系列 3 = 民間消費
%   系列 4 = 輸出
%   系列 5 = 為替レート

%  sign_restrictionの設定
%  Imp_t1(系列数)

%if Imp_t(1,1:1)>=0 & Imp_t(2,1:1)>=0 & Imp_t(5,1:1)>= 0 ;
   if Imp_t(1,1:1)>=0 & Imp_t(2,1:1)>=0 & Imp_t(5,1:1)>= 0 ; 
    
    % 財政政策ショックに対して政府消費(系列数;1)の反応は 1期から3期 までプラス
    % 財政政策ショックに対してGDP(系列数;2)の反応は 1期から3期 までプラス
    % 財政政策ショックに対して為替レート(系列数;5)の反応は 1期から3期 までプラス
    
    d=1;
else
    d=0;
end
end
end