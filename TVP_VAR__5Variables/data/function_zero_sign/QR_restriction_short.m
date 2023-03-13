% function [P] = QR_restriction_short(Sigma)

%　Rubio-Ramirez, Waggoner, and Zha (2010, RES, p665-696)
%  
%  p.678  Sec.5.2 A monetary SVAR

% 制約条件の設定
%  0の制約数　4
Q1 = [0 1 0 0 0; ...
      0 0 1 0 0; ...
      0 0 0 1 0; ...
      0 0 0 0 1;
      0 0 0 0 0];

%  0の制約数　3
Q2 = [0 0 1 0 0; ...
      0 0 0 1 0; ...
      0 0 0 0 1; ...
      0 0 0 0 0 ];
  
% 0の 制約数　3 (しかし2個が望ましい )  
Q3 = [1 0 0 0 0 ; ...
      0 1 0 0 0 ; ...
      0 0 0 0 1 ; ...
      0 0 0 0 0 ];

% 0の 制約数　1  
Q4 = [0 0 0 0 1; ...
      0 0 0 0 0; ...
      0 0 0 0 0; ...
      0 0 0 0 0 ];
  
% 0の 制約数　0    
Q5 = [0 0 0 0 0; ...
      0 0 0 0 0; ...
      0 0 0 0 0; ...
      0 0 0 0 0 ];


% 誘導系VARモデルの分散共分散行列
Sigma = [1   0.5  0.5  0.5 0.2;...
        0.5 4.25 2.5  0.5  0.25;...
        0.5  2.5  3   0.5  0.1;...
        0.5  0.5  0.5  2   0.15; ...
        0.2  0.25  0.1 0.15   4 ] ;

 invsig = inv(Sigma);
 A0 = chol(invsig);

 % inv(A0'*A0);
% invA0 = inv(A0');
% invA0'*invA0*A0'*A0; %
% invA0'*invA0;
% chol(sigma)'*chol(sigma);

%---------------------------
%disp('短期制約');
F = [A0];

% 行が0だけのものを削除
Q1_bar = abs(orth(Q1'));
Q2_bar = abs(orth(Q2'));
Q3_bar = abs(orth(Q3'));
Q4_bar = abs(orth(Q4'));

% 制約 1 の処理 
Q1_til = (Q1_bar' * F);
% QR decomposition
[Q,T]=qr(Q1_til');
P1 = Q(:,5)'; % 直交行列 Pの1列目

% 制約 2 の処理 
Q2_til = [ Q2_bar'*F; P1];
% QR decomposition
[Q,T]=qr(Q2_til');
P2 = Q(:,5)';  % 直交行列 Pの2列目

% 制約 3 の処理 
Q3_til = [ Q3_bar'*F; P1; P2];
% QR decomposition
[Q,T]=qr(Q3_til');
P3 = Q(:,5)';  % 直交行列 Pの2列目

% 制約 4 の処理 
Q4_til = [ Q4_bar'*F; P1; P2; P3];
% QR decomposition
[Q,T]=qr(Q4_til');
P4 = Q(:,5)';  % 直交行列 Pの2列目

% 制約 5 の処理
Q5_til = [P1; P2; P3; P4] ;
% QR decomposition
[Q,T]=qr(Q5_til');
P5 = Q(:,5)'; % 直交行列 Pの3列目

%disp('直交行列(ユニタリ行列) i.e.,　P*P = I');
P=[P1' P2' P3' P4' P5'];  % 直交行列、ユニタリ行列　P'*P = I

%disp('直交行列の検算 P*P = I');
disp(' ')
disp('P^T * P = ')
P'*P  %　直交行列の検算

% 制約の検算
disp( ' ' )
disp('A0 = Sig * P = ')
A0*P
