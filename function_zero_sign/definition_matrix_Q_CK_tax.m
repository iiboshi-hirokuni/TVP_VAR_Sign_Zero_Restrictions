
% 制約条件の設定
% 短期と長期のインパルス応答の制約 1 --> 1の数は4つ (インパルス応答が0となる制約)
%empty
Q1  =  [0 0 0 0 0; ...   % short run y
       0 0 0 0 0; ...   % tax
       0 0 0 0 0;...    % g
       0 0 0 0 0; ...  % int
       0 0 0 0 0; ...  %  pi
      0 0 0 0 0; ...   % short y    
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0; ... % pi  
      0 0 0 0 0; ...   % long run  y 
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0]' ;  % pi
  
% 短期と長期のインパルス応答の制約 2  --> 1の数は3つ
%empty
Q2  =  [0 0 0 0 0; ...   % short run y
       0 0 0 0 0; ...   % tax
       0 0 0 0 0;...    % g
       0 0 0 0 0; ...  % int
       0 0 0 0 0; ...  %  pi
      0 0 0 0 0; ...   % short y    
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0; ... % pi  
      0 0 0 0 0; ...   % long run  y 
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0]' ;  % pi

 % 短期と長期のインパルス応答の制約 3 --> 1の数は2つ
 % Tax adjusting Rule
Q3 = [0 0 0 0 0; ...   % short run y
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0; ...  %  pi
      0 0 0 0 0; ...   % short y 
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0; ... % pi  
      1 0 0 0 0; ...   % long run  y
      0 1 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0 ]' ;  % pi
 
  % 短期と長期のインパルス応答の制約 4 --> 1の数は1つ
  %empty
  Q4 =  [0 0 0 0 0; ...   % short run y
       0 0 0 0 0; ...   % tax
       0 0 0 0 0;...    % g
       0 0 0 0 0; ...  % int
       0 0 0 0 0; ...  %  pi
      0 0 0 0 0; ...   % short y    
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0; ... % pi  
      0 0 0 0 0; ...   % long run  y 
      0 0 0 0 0; ...   % tax
      0 0 0 0 0;...    % g
      0 0 0 0 0; ...  % int
      0 0 0 0 0]' ;  % pi


  Q=zeros(nk,3*nk,4);
Q(:,:,1) = Q1;
Q(:,:,2) = Q2;
Q(:,:,3) = Q3;
Q(:,:,4) = Q4;
  
