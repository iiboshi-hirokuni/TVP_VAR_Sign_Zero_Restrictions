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

function [mimp, sign_t,A0_t] = impulse_zero_sign_v2...
    (nl, nlen, mb, ma, mh, method_type, policy_type,svar_type)

% global method_type policy_type ;

% disp(policy_type)

ns = size(mh, 1);
nk = size(mh, 2);

% mimp = zeros(ns*nlen, nk^2);
mimp = zeros(nlen, nk^2,ns);

amOmsq = zeros(nk, nk, ns);
my = zeros(nl+nlen, nk);
mbs = [mb ; (ones(nlen, 1) * mb(ns, :))];

save_P=zeros(nk,nk,ns);

sign_P = ones(nk,ns);
sign_chk = zeros(ns,1);
sign_t = zeros(ns,1);
A0_t   = zeros(ns,nk);

vh = mean(mh(nl+1:end, :));     % average shock size
% save('tvpvar_imp_vh.xls', 'vh',  '-ascii', '-tabs');  % save average shock size

 for i = 1:nk 

  for t = nl+1 : ns  
    if i == 1 
       amOmsq(:, :, t) = real((fAt(ma(t, :), nk))\diag(exp(vh/2)));
%        amOmsq(:, :, t) = inv(fAt(ma(t, :), nk)) * diag(exp(vh/2));
    
        switch method_type
         case 'Arias_et_al'
           switch policy_type
               case 'bench_mark'                   
                   [save_P(:,:,t),sign_P(:,t),sign_chk(t),A0_t(t,:)] = ...
                       zero_sign_restriction_v4(amOmsq,nk, t, nl,nlen, mbs, svar_type);
               case 'Choleski'                     
                   save_P(:,:,t) = eye(nk);
                   sign_P(:,t) =1;
                   sign_chk(t)=1;
           end
         
        end   
    end

 
  if (sign_chk(t) == 1) 
      P = squeeze(save_P(:,:,t));
         Phi0= reshape(mbs(t,:)',nk*nl,nk);
     %%   
%          for h=1:nlen-1
%              my(h,:) = impulseresponce(Phi0,amOmsq(:,:,t),sign_P(i,t)*P(:,i),h);   
%          end 
%             mimp(:,(i-1)*nk+1:i*nk,t)= my(1:nlen, :);
      %%     
        my(nl+1, :) = (  sign_P(i,t)*amOmsq(:,:,t)*P(:,i) )' ;         
        for j = nl+2 : nl+nlen
            my(j, :) = mbs(t+j-nl-1,:) * fXt(my(j-nl:j-1,:), 0)';
        end      
        mimp(:,(i-1)*nk+1:i*nk,t)= my(nl+1:end, :);
     %%
       sign_t(t) = 1;
     
  else % (sign_chk(t) == 0) 
     mimp(:,(i-1)*nk+1:i*nk,t)= zeros(nlen, nk);
     sign_t(t) = 0;
  end
  
   
  end
  
end


		
