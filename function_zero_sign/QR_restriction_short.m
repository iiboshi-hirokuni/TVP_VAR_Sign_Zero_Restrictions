% function [P] = QR_restriction_short(Sigma)

%�@Rubio-Ramirez, Waggoner, and Zha (2010, RES, p665-696)
%  
%  p.678  Sec.5.2 A monetary SVAR

% ��������̐ݒ�
%  0�̐��񐔁@4
Q1 = [0 1 0 0 0; ...
      0 0 1 0 0; ...
      0 0 0 1 0; ...
      0 0 0 0 1;
      0 0 0 0 0];

%  0�̐��񐔁@3
Q2 = [0 0 1 0 0; ...
      0 0 0 1 0; ...
      0 0 0 0 1; ...
      0 0 0 0 0 ];
  
% 0�� ���񐔁@3 (������2���]�܂��� )  
Q3 = [1 0 0 0 0 ; ...
      0 1 0 0 0 ; ...
      0 0 0 0 1 ; ...
      0 0 0 0 0 ];

% 0�� ���񐔁@1  
Q4 = [0 0 0 0 1; ...
      0 0 0 0 0; ...
      0 0 0 0 0; ...
      0 0 0 0 0 ];
  
% 0�� ���񐔁@0    
Q5 = [0 0 0 0 0; ...
      0 0 0 0 0; ...
      0 0 0 0 0; ...
      0 0 0 0 0 ];


% �U���nVAR���f���̕��U�����U�s��
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
%disp('�Z������');
F = [A0];

% �s��0�����̂��̂��폜
Q1_bar = abs(orth(Q1'));
Q2_bar = abs(orth(Q2'));
Q3_bar = abs(orth(Q3'));
Q4_bar = abs(orth(Q4'));

% ���� 1 �̏��� 
Q1_til = (Q1_bar' * F);
% QR decomposition
[Q,T]=qr(Q1_til');
P1 = Q(:,5)'; % �����s�� P��1���

% ���� 2 �̏��� 
Q2_til = [ Q2_bar'*F; P1];
% QR decomposition
[Q,T]=qr(Q2_til');
P2 = Q(:,5)';  % �����s�� P��2���

% ���� 3 �̏��� 
Q3_til = [ Q3_bar'*F; P1; P2];
% QR decomposition
[Q,T]=qr(Q3_til');
P3 = Q(:,5)';  % �����s�� P��2���

% ���� 4 �̏��� 
Q4_til = [ Q4_bar'*F; P1; P2; P3];
% QR decomposition
[Q,T]=qr(Q4_til');
P4 = Q(:,5)';  % �����s�� P��2���

% ���� 5 �̏���
Q5_til = [P1; P2; P3; P4] ;
% QR decomposition
[Q,T]=qr(Q5_til');
P5 = Q(:,5)'; % �����s�� P��3���

%disp('�����s��(���j�^���s��) i.e.,�@P*P = I');
P=[P1' P2' P3' P4' P5'];  % �����s��A���j�^���s��@P'*P = I

%disp('�����s��̌��Z P*P = I');
disp(' ')
disp('P^T * P = ')
P'*P  %�@�����s��̌��Z

% ����̌��Z
disp( ' ' )
disp('A0 = Sig * P = ')
A0*P
