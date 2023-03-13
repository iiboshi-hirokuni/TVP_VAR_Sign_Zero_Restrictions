%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  mimp = impulse(nl, nlen, mb, ma, mh)
%%
%%  "impulse" computes time-varying impulse response
%%  (using average shock size)
%%
%%  [input]
%%    nl:       # of lags
%%    nlen:     length of response to compute
%%    mb,ma,mh: time-varying parameters
%%
%%  [output]
%%    mimp:  (ns*nlen)*(nk^2) matrix
%%

function mimp = impulse(nl, nlen, mb, ma, mh)

% sign-restriction　のスイッチ 
sign_res = 0;   %  スイッチ　オン=1,  オフ=0 
zero_res = 1;   %  スイッチ　オン=1,  オフ=0 　% Dec/23/2014 追加

ns = size(mh, 1);
nk = size(mh, 2);

% n = 200; % # of draw for sign restriction       Aug/3/2013

mimp = zeros(ns*nlen, nk^2);
amOmsq = zeros(nk, nk, ns);
my = zeros(nl+nlen, nk);
mbs = [mb ; (ones(nlen, 1) * mb(ns, :))];

save_P=zeros(nk,nk,ns);

vh = mean(mh(nl+1:end, :));     % average shock size
save('tvpvar_imp_vh.xls', 'vh',  '-ascii', '-tabs');  % save average shock size

 for i = 1:nk 
%  for i = 1:4:5  % Nov 30, 2013 for MFdiagonalization

  for t = nl+1 : ns  
    if i == 1 
       amOmsq(:, :, t) = inv(fAt(ma(t, :), nk)) ...
                       * diag(exp(vh/2));
       if (sign_res == 1)   %　挿入 Dec/10/2013
             [save_P(:,:,t),FP(t),MP(t),BC(t)] = sign_restriction(amOmsq,nk, t, nl,nlen,mbs);  %　挿入 Dec/10/2013
       elseif(zero_res == 1)
             [save_P(:,:,t),FP(t),MP(t),BC(t)] = zero_sign_restriction(amOmsq,nk, t, nl,nlen,mbs);  %　挿入 Dec/23/2014   
             
       end           
    end 

 %  start  sign_restriction
 
 % S = zeros(nk,nk,n);      % Aug/4/2013
 % P = zeros(nk,nk);        % Aug/4/2013
 % T = zeros(nk,nk);        % Aug/4/2013
 
 if (sign_res == 1)
     if(i==1);        % sign_restriction ON        %　挿入 July/18/2013/
     
  % for k = 1 : n           % Aug/3/2013 
  %   S(:,:,k) = sign_restriction(amOmsq,nk, t, nl,nlen,mbs); % Aug/4/2013
  % end                     % Aug/4/2013

    % T = sum(S,3);         % Aug/4/2013
    % P = (1/n).*T          % Aug/4/2013

        P = squeeze(save_P(:,:,t));
         my(nl+1, :) = (  FP(t)*amOmsq(:,:,t)*P(:,i) )' ;           %　挿入 July/18/2013/

     elseif (i==6);  % Dec/10/2013 
     
          P = squeeze(save_P(:,:,t));          
          my(nl+1, :) = (  MP(t)*amOmsq(:,:,t)*P(:,i) )' ;  
    
     elseif (i==2);  % Dec/13/2013 
     
          P = squeeze(save_P(:,:,t));          
           my(nl+1, :) = (  BC(t)*amOmsq(:,:,t)*P(:,i) )' ;  
     
     else    % Dec/10/2013 
          P = squeeze(save_P(:,:,t));
           my(nl+1, :) = (  FP(t)*amOmsq(:,:,t)*P(:,i) )' ;       
     end
     
     
 elseif (zero_res == 1)   % Dec/23/2014 追加
     if(i==1);        % Zero_restriction ON    
        P = squeeze(save_P(:,:,t));
         my(nl+1, :) = (  FP(t)*amOmsq(:,:,t)*P(:,i) )' ;           
     elseif (i==3);        
          P = squeeze(save_P(:,:,t));          
          my(nl+1, :) = (  MP(t)*amOmsq(:,:,t)*P(:,i) )' ;  
     elseif (i==2);        
          P = squeeze(save_P(:,:,t));          
          my(nl+1, :) = (  BC(t)*amOmsq(:,:,t)*P(:,i) )' ;       
     else    
          P = squeeze(save_P(:,:,t));
         my(nl+1, :) = (  FP(t)*amOmsq(:,:,t)*P(:,i) )' ;     
     end           
     
 else %if  sign_res == 0;     %     sign_restriction OFF 
          my(nl+1, :) = amOmsq(:, i, t)'; 
         
 end  
    
    for j = nl+2 : nl+nlen
      my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
    end
    
    mimp((t-1)*nlen+1:t*nlen, (i-1)*nk+1:i*nk) ...
      = my(nl+1:end, :);
  end  
  
  
end


		
