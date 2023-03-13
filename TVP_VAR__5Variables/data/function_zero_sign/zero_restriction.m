function [P, FP, MP, BC] = zero_restriction(amOmsq,nk, t, nl,nlen,mbs)

%  2014/Dec/23 �쐬

%�@Rubio-Ramirez, Waggoner, and Zha (2010) RES, p665-696
%  
%  p.679  Sec.5.3. Impulse Response 

% ��������̐ݒ�
% �Z���ƒ����̃C���p���X�����̐��� 1 --> 1�̐���5�� (�C���p���X������0�ƂȂ鐧��)
% ��������V���b�N
Q1 = [0 0 0 0 0 0; ...   % short run g
      0 0 0 0 0 0; ...   % y
      0 0 0 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      0 0 0 0 0 0; ...  % pi
      0 0 0 0 1 0; ...  % int
      1 0 0 0 0 0; ...   % long run  g
      0 1 0 0 0 0; ...   % y
      0 0 1 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      0 0 0 1 0 0; ...  % pi
      0 0 0 0 0 0  ]' ;     % int
  
% �Z���ƒ����̃C���p���X�����̐��� 2  --> 1�̐���4��
% �i�C�z��
Q2 =  [0 0 0 0 0 0; ...   % short run g
       0 0 0 0 0 0; ...   % y
       0 0 0 0 0 0; ...  % cons
       0 0 0 0 0 0;...     % dbt
       0 0 0 0 0 0; ...  % pi
       0 0 0 0 0 0; ...  % int
       0 0 0 0 0 0; ...   % long run  g
       1 0 0 0 0 0; ...   % y
       0 1 0 0 0 0; ...  % cons
       0 0 0 0 0 0;...     % dbt
       0 0 1 0 0 0; ...  % pi
       0 0 0 1 0 0  ]' ;     % int

 % �Z���ƒ����̃C���p���X�����̐��� 3 --> 1�̐���3��
 % ���Z���� 
Q3 = [0 0 1 0 0 0; ...   % short run g
      0 0 0 0 0 0; ...   % y
      0 0 0 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      0 0 0 0 0 0; ...  % pi
      0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0; ...   % long run  g
      0 0 0 0 0 0; ...   % y
      0 0 0 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      0 1 0 0 0 0; ...  % pi
      1 0 0 0 0 0  ]' ;     % int
 
  % �Z���ƒ����̃C���p���X�����̐��� 4 --> 1�̐���2��
  %�@�����V���b�N
Q4 =  [0 0 0 0 0 0; ...   % short run g
       0 0 0 0 0 0; ...   % y
      0 0 0 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      1 0 0 0 0 0; ...  % pi
      0 1 0 0 0 0; ...  % int
      0 0 0 0 0 0; ...   % long run  g
      0 0 0 0 0 0; ...   % y
      0 0 0 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      0 0 0 0 0 0; ...  % pi
      0 0 0 0 0 0  ]' ;     % intt

 % �Z���ƒ����̃C���p���X�����̐��� 5 --> 1�̐���1��
 % ���v�V���b�N
Q5  =  [0 0 0 0 0 0; ...   % short run g
       0 0 0 0 0 0; ...   % y
      0 0 0 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      0 0 0 0 0 0; ...  % pi
      0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0; ...   % long run  g
      1 0 0 0 0 0; ...   % y
      0 0 0 0 0 0; ...  % cons
      0 0 0 0 0 0;...     % dbt
      0 0 0 0 0 0; ...  % pi
      0 0 0 0 0 0  ]' ;     % int

% �U���nVAR���f���̕��U�����U�s��
%sigma = [1   0.5  1  ;...
%         0.5 4.25 2.5;...
%         1   2.5  3  ];
% invsig = inv(sigma);
% A0 = chol(invsig);
% inv(A0'*A0);
% invA0 = inv(A0');
% invA0'*invA0*A0'*A0; %
% invA0'*invA0;
% chol(sigma)'*chol(sigma);

%---------------------------

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

% �s��0�����̂��̂��폜
Q1_bar = abs(orth(Q1'));
Q2_bar = abs(orth(Q2'));
Q3_bar = abs(orth(Q3'));
Q4_bar = abs(orth(Q4'));
Q5_bar = abs(orth(Q5'));

% ���� 1 �̏��� 
Q1_til = (Q1_bar' * F);
% QR decomposition
[Q,T]=qr(Q1_til');
P1 = Q(:,6)'; % �����s�� P��1���

% ���� 2 �̏��� 
Q2_til = [ Q2_bar'*F; P1];
% QR decomposition
[Q,T]=qr(Q2_til');
P2 = Q(:,6)';  % �����s�� P��2���

% ���� 3 �̏��� 
Q3_til = [ Q3_bar'*F; P1; P2];
% QR decomposition
[Q,T]=qr(Q3_til');
P3 = Q(:,6)';  % �����s�� P��2���

% ���� 4 �̏��� 
Q4_til = [ Q4_bar'*F; P1; P2; P3];
% QR decomposition
[Q,T]=qr(Q4_til');
P4 = Q(:,6)';  % �����s�� P��2���

% ���� 5 �̏��� 
Q5_til = [ Q5_bar'*F; P1; P2; P3; P4];
% QR decomposition
[Q,T]=qr(Q5_til');
P5 = Q(:,6)';  % �����s�� P��2���

% ���� 6 �̏���
Q6_til = [P1; P2; P3; P4; P5] ;
% QR decomposition
[Q,T]=qr(Q6_til');
P6 = Q(:,6)'; % �����s�� P��3���

%disp('�����s��(���j�^���s��) i.e.,�@P*P = I');
P=[P1' P2' P3' P4' P5' P6'];  % �����s��A���j�^���s��@P'*P = I

i = 1;    % ��������(���{�V���b�N)�V���b�N���N�����W����order
 my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;
 Imp_FP =  my(nl+1, :);
 
if( Imp_FP(1)>0)  % Gov >0
    FP = 1;      
elseif (-1*Imp_FP(1)>0 )
    FP = -1;     
end    

i = 2;    % �i�C�z�V���b�N���N�����W����order 
  my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;
 Imp_BC =  my(nl+1, :);
 
if( Imp_BC(2)>0)  % Y >0
    BC = 1;      
elseif (-1*Imp_BC(2)>0 )
    BC = -1;     
end    

 i = 3;    % ���Z����V���b�N���N�����W����order 
   my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;
 Imp_MP =  my(nl+1, :);
 
if( Imp_MP(6)>0)   % r >0
    MP = 1;      
elseif (-1*Imp_MP(6)>0 )
    MP = -1;     
end   

 


%========================================================== 
%  �[������̎����̎��{
% =========================================================

test = 0;  % yes--> 1 no--> 0

if test == 1

if t== 15 
disp('�����s��̌��Z P^T * P = I');
 disp( P'*P );  %�@�����s��̌��Z�@�P�ʍs��ɂȂ邩
% 
% % ����P�ƃV���b�N1�̌��Z 
 e1 = [1; 0; 0; 0; 0; 0 ]; % �V���b�N1
disp('����P�ƃV���b�N1�̌��Z: Q1*F*P*e1, �����s�� P ���� ');
%   disp( Q1*F*P*e1 );  % �����s�� P ����@���ׂĂ̒l��0�ɋ߂��Ȃ邩
  disp( Q1*F*P(:,1) ); 
% % ����2�ƃV���b�N2�̌��Z 
 e2 = [0;1;0; 0; 0; 0]; % �V���b�N2
disp('����2�ƃV���b�N2�̌��Z: Q2*F*P*e2, �����s�� P ���� ');
% disp( Q2*F*P*e2 );  % �����s�� P ����
 disp( Q2*F*P(:,2) );  
% % ����3�ƃV���b�N3�̌��Z 
 e3 = [0;0;1; 0; 0; 0]; % �V���b�N2
disp('����3�ƃV���b�N3�̌��Z: Q2*F*P*e3, �����s�� P ���� ');
% disp( Q3*F*P*e3 );   % �����s�� P ����
  disp( Q3*F*P(:,3) );  
% % ����4�ƃV���b�N4�̌��Z 
 e4 = [0;0;0; 1; 0; 0]; % �V���b�N2
disp('����4�ƃV���b�N4�̌��Z: Q4*F*P*e4, �����s�� P ���� ');
% disp( Q4*F*P*e4 );   % �����s�� P ����
   disp( Q4*F*P(:,4) ); 
% % ����5�ƃV���b�N5�̌��Z 
 e5 = [0;0;0; 0; 1; 0]; % �V���b�N2
disp('����5�ƃV���b�N5�̌��Z: Q2*F*P*e5, �����s�� P ���� ');
% disp( Q5*F*P*e5 );   % �����s�� P ����
 disp( Q5*F*P(:,5) );   
end

end

end