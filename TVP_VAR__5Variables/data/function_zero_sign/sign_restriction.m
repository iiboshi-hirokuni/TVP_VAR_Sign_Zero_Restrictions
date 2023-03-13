function [P,FP,MP,BC,chk] = sign_restriction(amOmsq,nk, t, nl,nlen,mbs)

chk =1;

% nk = �n��̐�
% amOmsq = inv(A)*Sigma = 1���ڂ̃��X�|���X����
% t = �V���b�N��������period�@(nl+1(���O�̎���+1)�`ns(�f�[�^�̃T���v������)) 

a=0; % Nov30, 2013
b=0; % Nov30, 2013
d=0;

my = zeros(nl+nlen, nk);   % �C���p���X���X�|���X

while d==0

% Omega = QR decompotion
% RWZ (2010) p688, algorithm 2 

x = randn(nk,nk);

[P, R]=qr(x);        %  i.e., P * P' = I

i = 1;    % ��������(���{�V���b�N)�V���b�N���N�����W����order 
%  my(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1���ڂ̐���t�����X�|���X
my(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;

%  2���ȍ~�́@�C���p���X����
  for j = nl+2 : nl+nlen
       my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
  end
  
  Imp_FP = zeros(nk,nlen);
  Imp_FP =  my(nl+1: nl+nlen,:)' ;

%  �n�� 1 = gov cons & inv
%  �n�� 2 = real GDP
%  �n�� 3 = private consumption
%  �n�� 4 = Debt GDP ratio
%  �n�� 5 = GDP deflator
%  �n�� 6 = Interest Rate
%  �n�� 7 = real exchange rate
%  �n�� 8 = net export

%  sign_restriction�̐ݒ�
%  Imp_t1(�n��)

if( Imp_FP(1,1:1)>0) &(Imp_FP(2,1:1)>0) &(Imp_FP(4,1:1)>0);
    FP = 1;     a=1; % Nov30, 2013
elseif (-1*Imp_FP(1,1:1)>0 )&(-1*Imp_FP(2,1:1)>0)&(-1*Imp_FP(4,1:1)>0); % Dec 11, 2013
    FP = -1;    a=1; % Dec 11, 2013

    
    % ��������V���b�N�ɑ΂��Đ��{�x�o(�n��;1)�̔����� 1���v���X
    % ��������V���b�N�ɑ΂���GDP(�n��;2)�̔����� 1���v���X
    % ��������V���b�N�ɑ΂���DEBT(�n��;4)�̔����� 1���v���X   
    
else
    a=0; % Nov30, 2013
end

% Nov30, 2013 %%%%%%%%%%

% a=1; FP = 1;

if a==0
    d=0;
    
else     
    
i = 3;    % ���Z����(���q���V���b�N)�V���b�N���N�����W����order 
%  my2(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1���ڂ̐���t�����X�|���X
 my2(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;

%  2���ȍ~�́@�C���p���X����
  for j = nl+2 : nl+nlen
       my2(j, :) = mbs(t+j-nl-1,:) * fXt(my2(j-nl:j-1,:), 0)';
  end
  
  Imp_MP = zeros(nk,nlen);
  Imp_MP =  my2(nl+1: nl+nlen,:)' ;

    
 if (Imp_MP(2,1:1)<0 )&( Imp_MP(5,1:1)<0 )&( Imp_MP(6,1:1)>0 );
                     MP = 1;      b=1; 
elseif (-1*Imp_MP(2,1:1)<0)& (-1*Imp_MP(5,1:1)<0)& (-1*Imp_MP(6,1:1)>0);
      MP = -1;     b=1;    

    % ���Z�����߃V���b�N�ɑ΂���GDP(�n��;2)�̔����� 1���}�C�i�X
    % ���Z�����߃V���b�N�ɑ΂��ăC���t����(�n��;5)�̔����� 1���}�C�i�X
    % ���Z�����߃V���b�N�ɑ΂��ė��q��(�n��;6)�̔����� 1���v���X
  
else
    b=0;
 end
 
% Dec13, 2013 %%%%%%%%%%

% a=1; FP = 1;

if b==0
    d=0;
    
else     
    
i = 2;    % ���Y���V���b�N���N�����W����order 
%  my3(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1���ڂ̐���t�����X�|���X
 my3(nl+1, :) = (  amOmsq(:,:,t)*P(:,i) )' ;

%  2���ȍ~�́@�C���p���X����
  for j = nl+2 : nl+nlen
       my3(j, :) = mbs(t+j-nl-1,:) * fXt(my3(j-nl:j-1,:), 0)';
  end
  
  Imp_BC = zeros(nk,nlen);
  Imp_BC =  my3(nl+1: nl+nlen,:)' ;

if (Imp_BC(2,1:1)>0 )&( Imp_BC(3,1:1)>0 )&( Imp_BC(4,1:1)<0 );
                     BC = 1;      d=1; 
elseif (-1*Imp_BC(2,1:1)>0)&(-1*Imp_BC(3,1:1)>0 )& (-1*Imp_BC(4,1:1)<0);
      BC = -1;     d=1;    
    % ���Y���V���b�N�ɑ΂���GDP(�n��;2)�̔����� 1���v���X
    % ���Y���V���b�N�ɑ΂��ď���(�n��;3)�̔����� 1���v���X
    % ���Y���V���b�N�ɑ΂���DEBT(�n��;4)�̔����� 1���}�C�i�X
  
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

