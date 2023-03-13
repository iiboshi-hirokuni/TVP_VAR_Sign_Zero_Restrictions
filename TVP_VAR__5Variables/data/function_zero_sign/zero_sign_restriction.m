function [P, FP, MP, BC, sign_chk] = zero_sign_restriction(amOmsq,nk, t, nl,nlen,mbs)

%�@Rubio-Ramirez, Waggoner, and Zha (2010) RES, p665-696
%  p.679  Sec.5.3. Impulse Response 

% ��������̐ݒ�
% �Z���ƒ����̃C���p���X�����̐��� 1 --> 1�̐���6�� (�C���p���X������0�ƂȂ鐧��)
% ��������V���b�N
Q1 = [0 0 0 0 0 0 0; ...   % short run y
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0; ...  %  pi
      0 0 0 0 0 0 0; ...   % long run  y    
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0 ]' ;  % pi
  
% �Z���ƒ����̃C���p���X�����̐��� 2  --> 1�̐���5��
% �i�C�z��
Q2 =  [0 0 0 0 0 0 0; ...   % short run y
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0; ...  %  pi
      0 0 0 0 0 0 0; ...   % long run  y    
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0 ]' ;  % pi

 % �Z���ƒ����̃C���p���X�����̐��� 3 --> 1�̐���4��
 % ���Z���� 
Q3 = [0 0 0 0 0 0 0; ...   % short run y
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0; ...  %  pi
      0 0 0 0 0 0 0; ...   % long run  y    
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0 ]' ;  % pi
 
  % �Z���ƒ����̃C���p���X�����̐��� 4 --> 1�̐���3��
  %�@�����V���b�N
Q4 =  [0 0 0 0 0 0 0; ...   % short run y
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0; ...  %  pi
      0 0 0 0 0 0 0; ...   % long run  y    
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0 ]' ;  % pi

 % �Z���ƒ����̃C���p���X�����̐��� 5 --> 1�̐���2��
 % ���v�V���b�N
Q5  =  [0 0 0 0 0 0 0; ...   % short run y
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0; ...  %  pi
      0 0 0 0 0 0 0; ...   % long run  y    
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0 ]' ;  % pi
  
 % �Z���ƒ����̃C���p���X�����̐��� 6 --> 1�̐���1��
 % ���v�V���b�N
Q6  =  [0 0 0 0 0 0 0; ...   % short run y
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0; ...  %  pi
      0 0 0 0 0 0 0; ...   % long run  y    
      0 0 0 0 0 0 0; ...   % c
      0 0 0 0 0 0 0; ...   % tax
      0 0 0 0 0 0 0;...    % g
      0 0 0 0 0 0 0; ...  % inv
      0 0 0 0 0 0 0; ...  % int
      0 0 0 0 0 0 0 ]' ;  % pi 
  
  

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

%====================================================
%  Zero Restriction
%====================================================

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
Q6_bar = abs(orth(Q6'));

% ���� 1 �̏��� 
Q1_til = (Q1_bar' * F);
% QR decomposition
[Q,T]=qr(Q1_til');
P1 = Q(:,7)'; % �����s�� P��1���

% ���� 2 �̏��� 
Q2_til = [ Q2_bar'*F; P1];
% QR decomposition
[Q,T]=qr(Q2_til');
P2 = Q(:,7)';  % �����s�� P��2���

% ���� 3 �̏��� 
Q3_til = [ Q3_bar'*F; P1; P2];
% QR decomposition
[Q,T]=qr(Q3_til');
P3 = Q(:,7)';  % �����s�� P��2���

% ���� 4 �̏��� 
Q4_til = [ Q4_bar'*F; P1; P2; P3];
% QR decomposition
[Q,T]=qr(Q4_til');
P4 = Q(:,7)';  % �����s�� P��2���

% ���� 5 �̏��� 
Q5_til = [ Q5_bar'*F; P1; P2; P3; P4];
% QR decomposition
[Q,T]=qr(Q5_til');
P5 = Q(:,7)';  % �����s�� P��2���

% ���� 6 �̏��� 
Q6_til = [ Q5_bar'*F; P1; P2; P3; P4; P5];
% QR decomposition
[Q,T]=qr(Q5_til');
P6 = Q(:,7)';  % �����s�� P��2���

% ���� 6 �̏���
Q7_til = [P1; P2; P3; P4; P5; P6] ;
% QR decomposition
[Q,T]=qr(Q7_til');
P7 = Q(:,7)'; % �����s�� P��3���

%disp('�����s��(���j�^���s��) i.e.,�@P*P = I');
P=[P1' P2' P3' P4' P5' P6' P7'];  % �����s��A���j�^���s��@P'*P = I


%====================================================
%
%  Sign Restriction
%
%====================================================

sign_OK = 0;   BC=0;MP=0;  FP = 0;  Tax =0;


 i = 1;    % �i�C�z�V���b�N���N�����W����order 
   my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;
   Imp_BC =  my(nl+1, :);
 
 if( Imp_BC(1)>0) && (Imp_BC(2)>0) && (Imp_BC(3)>0) &&(Imp_BC(5)>0)  % Y >0, C>0, Tax>0, Inv>0
      BC = 1;      
 elseif (-1*Imp_BC(1)>0) &&(-Imp_BC(2)>0) &&( -Imp_BC(3)>0) && (-Imp_BC(5)>0 )
      BC = -1;  
 else 
     sign_OK = 0;  
 end    
end 
% 

if sign_OK == 1;
 i = 2;    % ���Z����V���b�N���N�����W����order 
   my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;
   Imp_MP =  my(nl+1, :);
 
 if( Imp_MP(6)>0)&&( Imp_MP(7)<0)   % r >0 & Pi < 0 
         MP = 1;      
  elseif (-1*Imp_MP(6)>0 )&&(-1*Imp_MP(7)<0)
        MP = -1;    
 else 
        sign_OK = 0;
  end
end  

if sign_OK == 1;
i = 3;    % ��������(���{�V���b�N)�V���b�N���N�����W����order
 my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;
 Imp_FP =  my(nl+1, :);

if( Imp_FP(4)>0)  % Gov >0 
    FP = 1;    sign_OK = 1;
elseif (-1*Imp_FP(4)>0 )
    FP = -1;    sign_OK = 1;  
end    

if sign_OK == 1;
i = 4;    % �ŋ��V���b�N���N�����W����order
 my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;
 Imp_Tax =  my(nl+1, :);

if( Imp_Tax(3)>0)  % tax >0
    Tax = 1;    sign_OK = 1;
elseif (-1*Imp_Tax(3)>0 )
    Tax = -1;    sign_OK = 1;  
end    


if sign_OK == 1;
     sign_chk = 1;
else
     sign_chk = 0;      
end

%========================================================== 
%  �[������̎����̎��{
% =========================================================
% 
% test = 0;  % yes--> 1 no--> 0
% 
% if test == 1
% 
% if t== 15 
% disp('�����s��̌��Z P^T * P = I');
%  disp( P'*P );  %�@�����s��̌��Z�@�P�ʍs��ɂȂ邩
% % 
% % % ����P�ƃV���b�N1�̌��Z 
%  e1 = [1; 0; 0; 0; 0; 0 ]; % �V���b�N1
% disp('����P�ƃV���b�N1�̌��Z: Q1*F*P*e1, �����s�� P ���� ');
% %   disp( Q1*F*P*e1 );  % �����s�� P ����@���ׂĂ̒l��0�ɋ߂��Ȃ邩
%   disp( Q1*F*P(:,1) ); 
% % % ����2�ƃV���b�N2�̌��Z 
%  e2 = [0;1;0; 0; 0; 0]; % �V���b�N2
% disp('����2�ƃV���b�N2�̌��Z: Q2*F*P*e2, �����s�� P ���� ');
% % disp( Q2*F*P*e2 );  % �����s�� P ����
%  disp( Q2*F*P(:,2) );  
% % % ����3�ƃV���b�N3�̌��Z 
%  e3 = [0;0;1; 0; 0; 0]; % �V���b�N2
% disp('����3�ƃV���b�N3�̌��Z: Q2*F*P*e3, �����s�� P ���� ');
% % disp( Q3*F*P*e3 );   % �����s�� P ����
%   disp( Q3*F*P(:,3) );  
% % % ����4�ƃV���b�N4�̌��Z 
%  e4 = [0;0;0; 1; 0; 0]; % �V���b�N2
% disp('����4�ƃV���b�N4�̌��Z: Q4*F*P*e4, �����s�� P ���� ');
% % disp( Q4*F*P*e4 );   % �����s�� P ����
%    disp( Q4*F*P(:,4) ); 
% % % ����5�ƃV���b�N5�̌��Z 
%  e5 = [0;0;0; 0; 1; 0]; % �V���b�N2
% disp('����5�ƃV���b�N5�̌��Z: Q2*F*P*e5, �����s�� P ���� ');
% % disp( Q5*F*P*e5 );   % �����s�� P ����
%  disp( Q5*F*P(:,5) );   
% end
% 
% end
% 
% end