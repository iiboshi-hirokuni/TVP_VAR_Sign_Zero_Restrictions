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

function [mimp, sign] = impulse_zero_sign(nl, nlen, mb, ma, mh)

ns = size(mh, 1);
nk = size(mh, 2);

% mimp = zeros(ns*nlen, nk^2);
mimp = zeros(nlen, nk^2,ns);

amOmsq = zeros(nk, nk, ns);
my = zeros(nl+nlen, nk);
mbs = [mb ; (ones(nlen, 1) * mb(ns, :))];

save_P=zeros(nk,nk,ns);

vh = mean(mh(nl+1:end, :));     % average shock size
save('tvpvar_imp_vh.xls', 'vh',  '-ascii', '-tabs');  % save average shock size

 for i = 1:nk 

  for t = nl+1 : ns  
    if i == 1 
       amOmsq(:, :, t) = inv(fAt(ma(t, :), nk)) ...
                       * diag(exp(vh/2));
                   
%       [save_P(:,:,t),FP(t),MP(t),BC(t),sign_chk(t)] = zero_sign_restriction(amOmsq,nk, t, nl,nlen,mbs);  %Å@ë}ì¸ Dec/23/2014   
             
%     [save_P(:,:,t),FP(t),MP(t),BC(t), sign_chk(t)] = sign_restriction(amOmsq,nk, t, nl,nlen,mbs);  %Å@
                   
%     [save_P(:,:,t),FP(t),MP(t),BC(t),sign_chk(t)] = zero_sign_restriction_1(amOmsq,nk, t, nl,nlen,mbs);  %Å@ë}ì¸ Dec/23/2014     
  
    [save_P(:,:,t),FP(t),MP(t),BC(t),sign_chk(t)] = sign_restriction_1(amOmsq,nk, t, nl,nlen,mbs);  %Å@ë}ì¸ Dec/23/2014     
   
    end          
   
  if (sign_chk(t) == 1) 
     if(i==1);        % ç‡ê≠ê≠çÙ    
        P = squeeze(save_P(:,:,t));
         my(nl+1, :) = (  FP(t)*amOmsq(:,:,t)*P(:,i) )' ;           
     elseif (i==3);       % ã‡óZê≠çÙ 
          P = squeeze(save_P(:,:,t));          
          my(nl+1, :) = (  MP(t)*amOmsq(:,:,t)*P(:,i) )' ;  
     elseif (i==2);        % åiãCèzä¬
          P = squeeze(save_P(:,:,t));          
          my(nl+1, :) = (  BC(t)*amOmsq(:,:,t)*P(:,i) )' ;       
     else    
          P = squeeze(save_P(:,:,t));
         my(nl+1, :) = (  FP(t)*amOmsq(:,:,t)*P(:,i) )' ;     
     end           
     
    
    for j = nl+2 : nl+nlen
      my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
    end    
    
     mimp(:,(i-1)*nk+1:i*nk,t)= my(nl+1:end, :);
     sign(t) = 1;
  else
     mimp(:,(i-1)*nk+1:i*nk,t)= zeros(nlen, nk);
     sign(t) = 0;
  end
  
 end  % 
    
end


		
