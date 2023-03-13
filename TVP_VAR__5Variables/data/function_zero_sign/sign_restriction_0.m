function [P] = sign_restriction(amOmsq,nk, t, nl,nlen,mbs)

% nk = �n��̐�
% amOmsq = inv(A)*Sigma = 1���ڂ̃��X�|���X����
% t = �V���b�N��������period�@(nl+1(���O�̎���+1)�`ns(�f�[�^�̃T���v������)) 

d=0;

my = zeros(nl+nlen, nk);   % �C���p���X���X�|���X

while d==0

%  Omega = QR decompotion
% RWZ (2010) p688, algorithm 2 

x = randn(nk,nk);

[P, R]=qr(x);        %  i.e., P * P' = I

i = 1;    % ��������(���{�V���b�N)�V���b�N���N�����W����order 
 my(nl+1, :) = ( P * amOmsq(:,i,t) )' ;  % 1���ڂ̐���t�����X�|���X

%  2���ȍ~�́@�C���p���X����
  for j = nl+2 : nl+nlen
       my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
  end
  
  Imp_t = zeros(nk,nlen);
  Imp_t =  my(nl+1: nl+nlen,:)' ;

% �@�n�� 1 = ���{����
%   �n�� 2 = real GDP
%   �n�� 3 = ���ԏ���
%   �n�� 4 = �A�o
%   �n�� 5 = �בփ��[�g

%  sign_restriction�̐ݒ�
%  Imp_t1(�n��)

%if Imp_t(1,1:1)>=0 & Imp_t(2,1:1)>=0 & Imp_t(5,1:1)>= 0 ;
   if Imp_t(1,1:1)>=0 & Imp_t(2,1:1)>=0 & Imp_t(5,1:1)>= 0 ; 
    
    % ��������V���b�N�ɑ΂��Đ��{����(�n��;1)�̔����� 1������3�� �܂Ńv���X
    % ��������V���b�N�ɑ΂���GDP(�n��;2)�̔����� 1������3�� �܂Ńv���X
    % ��������V���b�N�ɑ΂��Ĉבփ��[�g(�n��;5)�̔����� 1������3�� �܂Ńv���X
    
    d=1;
else
    d=0;
end
end
end