%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  malpha = ssmooth(my, amZ, amG2, mH2, va0, mH02)
%%
%%  "ssmooth" implements simulation smoother
%%  by de Jong & Shephard (1995)
%%
%%  [model]
%%      y_t         =  D_t + Z_t * alpha_t + v_t,   v_t ~ N(0, HH)              
%%      alpha_{t+1} =  T_t * alpha_t + R_t*e_t,     e_t ~ N(0, QQ)
%%      
%%      G_t*H_t' = O
%%
%%      y_t:     nk*1 vector
%%      Z_t:     nk*np matrix
%%      alpha_t: np*1 vector
%%      G_t:     nk*(nk+np) matrix
%%      H_t:     np*(nk+np) matrix
%%
%%  [input]
%%      my:     response (ns*nk vector)
%%      ma:
%%      mb:
%%      mh:
%%      amZ:    independent variable (nk*np*ns array)
%%      amG2:   G_t*G_t' (nk*nk*ns array)
%%      mH2:    H_t*H_t' = H*H' (np*np matrix)
%%      va0:    alpha_0 (np*1 vector)
%%      mH02:   H_0*H_0' (np*np matrix)
%%
%%  [output]
%%      malpha:  sampled state variable (np*ns matrix)
%%

 function [decomp_smooth amOmsq ]  = hist_decomp(my, ma, mb, mh, P, sign_P) %, amZ, amG2, mH2, va0, mH02)
% function [obsmean  amOmsq ] = hist_decomp(my, ma, mb, mh, P, sign_P)


%%--- set variables ---%%

ns = size(my, 1);    % # of time periods
nk = size(my, 2);    % # of series
% np = size(mH2, 1);    % # of state
nl = 4;              % # of lag
nstate = nl*nk;      % # of state

nshock = nk;         % # of shocks

nobs =ns;

vh = mean(mh(nl+1:end, :));     % average shock size
va = mean(ma(nl+1:end, :)); 

% va = va0;
% mP = mH02;
% vr = zeros(np, 1);
% mU = zeros(np);
% me = zeros(nk, ns);
% amDinv = zeros(nk, nk, ns);
% amL = zeros(np, np, ns);
% meta = zeros(np, ns);

obsmean   = zeros(nobs, nk);
% smooth1   = zeros(nobs, nk);
% smooth2   = zeros(nobs, nk);
obsvar    = zeros(nobs, nk);
nu_save   = zeros(nk,nobs);
% ft_save   = zeros(nk,nk*nobs); 
% kg_save   = zeros(nstate,nk*nobs);
at_save   = zeros(nobs, nstate);
% pt_save   = zeros(nstate, nstate*nobs);
% shock     = zeros(nobs, nshock);
st_save   = zeros(nobs,nk*nshock);

alpha = zeros(nk*nl,1);

ZZ = [eye(nk) zeros(nk,nk*(nl-1)) ];

   
DD = zeros(nk,1);
% HH = zeros(nk,nk);
HH = 0.01*eye(nk);
QQ = eye(nk);   % createcov(para(30:39,1));
VV = zeros(nshock,nk);

At = zeros(nstate,1); 
a1 = At;

% Pt = dlyap(TT,RR*QQ*RR');  
Pt = zeros(nstate);
% RR = [RR; zeros(nstate,nshock)];

% amOmsq = zeros(nk, nk, ns);
% mbs = [mb ; (ones(nlen, 1) * mb(ns, :))];
% compute likelihood with Kalman filter by Durbin & Koopman (2012, p85) 

t = 1;
while t <= nobs
    
%      amOmsq = inv(fAt(va, nk)) ; 
    amOmsq = inv(fAt(ma(t, :), nk))* diag(exp(mh(t, :)/2));
%                        * diag(exp(vh/2));
%                     * diag(exp(mh(t, :)/2));
                   
%    RR = ( sign_P(i)*amOmsq*P ) ;
    RR = ( amOmsq ) ;
%    RR = eye(nk);
   RR = [RR; zeros(nk*(nl-1),nk) ];
   
   
   TT = [ reshape( mb(t,:), nk:nl,nk)' ;...
          eye(nk) zeros(nk,nk*3); ...
          zeros(nk) eye(nk) zeros(nk, nk*2); ...
          zeros(nk, nk*2) eye(nk) zeros(nk, nk)];   

%  TT = [  0.1*eye(nk)  0.2*eye(nk)  0.205*eye(nk)  0.2*eye(nk); ...
%          eye(nk) zeros(nk,nk*3); ...
%          zeros(nk) eye(nk) zeros(nk, nk*2); ...
%          zeros(nk, nk*2) eye(nk) zeros(nk, nk)];   


   At1 = At;
   Pt1 = Pt;
   
   % Forecasting
   
   alphahat = TT*At1;
   Phat = TT*Pt1*TT' + RR*QQ*RR';
   yhat = ZZ*alphahat + DD;
   nu   = my(t,:) - yhat';
   
   Ft   = ZZ*Phat*ZZ' + HH + ZZ*RR*VV + (ZZ*RR*VV)';
   Ft   = 0.5*(Ft + Ft');
   
   K_g = TT*Phat*ZZ'*inv(Ft);   % Kalman Gain
  
   at_save(t,:) = alphahat;
   pt_save(:,(t-1)*nstate+1:t*nstate) = Phat;
   
%    loglh = loglh -0.5*size(my, 2)*log(2*pi)-0.5*log(det(Ft)) ...
%            - 0.5*nu*inv(Ft)*nu';
   
   % Updating
   At = alphahat + (Phat*ZZ' + RR*VV)*inv(Ft)*nu';
   Pt = Phat - (Phat*ZZ'+RR*VV)*inv(Ft)*(Phat*ZZ'+RR*VV)';
   
   %  store
%    obsmean(t,:) =yhat'  ; %alphahat(nstate/2,:)';
   
    obsmean(t,:) = ( ZZ*alphahat+DD)';
   
   obsvar(t,:)  = diag(Ft)';
   nu_save(:,t) = nu';
   Ft_save(:,(t-1)*nk+1:t*nk) = Ft;
   Kg_save(:,(t-1)*nk+1:t*nk) = K_g;
   
   t = t+1;
end  

% disturbance smoothing by Durbin & Koopman (2001, p??)
%                       by Durbin & Koopman (2012, p96)

r_t = zeros(nstate,1);
N_t = zeros(nstate,nstate); 
eta = zeros(nshock,nobs);

for t = nobs:-1:1 

  nu = nu_save(:,t);
  Ft =  Ft_save(:,(t-1)*nk+1:t*nk);
  K_g = Kg_save(:,(t-1)*nk+1:t*nk);
  Phat = pt_save(:,(t-1)*nstate+1:t*nstate);
  L_t = TT- K_g*ZZ; 
  eta(:,t) = QQ*RR'*r_t;
  
  r_t = ZZ'*inv(Ft)*nu + L_t'*r_t;     % Eq.(4.32)
  N_t = ZZ'*inv(Ft)*ZZ + L_t'*N_t*L_t; 
  
end

% /*   historical decomposition
% */

nvar    = nk;

for  index = 1:nshock;
%    st1 = zeros(40,1);
%    r1 = zeros(40,1);
%    r1(index) = r_t(index);  
   
%    st1 = a1 + p1*r_t;  %QQ(index,index)/sum(diag(QQ)) ;  % initialize    
    st1 = a1;                                             % initialize   

 for t =1:nobs-1;   
      st_save(t,1+(index-1)*nvar:index*nvar ) = (ZZ*st1)';
    
     shock_tmp = zeros(nshock,1);
     shock_tmp(index) = eta(index,t);     
     st = TT*st1 + RR*shock_tmp; 
     st1 = st;
 end;
 
end;

 decomp_smooth = st_save;

%====================================================================
% 
%   End of  historical decomposition
%
%==================================================================== 

% %--- Kalman filter ---%%
% 
% for i = 1 : ns
%     
%     TT = [reshape(mb(i,:),nk,nk*nl); ...
%           eye(nk) zeros(nk,nk*(nl-1));...
%           zeros(nk,nk)  eye(nk) zeros(nk,nk*(nl-2));...
%           zeros(nk,nk*2)  eye(nk) zeros(nk,nk) ];
%   
%    % forecasting   
%     alpha = TT*alpha;     
%     mD = ZZ* mP * ZZ' + mh(:,:,i);  % amG2
%     me(:, i) = my(i, :)' - (ZZ*alpha)';
%     
%     drD = rcond(mD);
%     if isnan(drD) || (drD < eps*10^2)
%         mDinv = eye(np) * 10;
%     else
%         mDinv = inv(mD);
%     end
%     amDinv(:, :, i) = mDinv;
% 
%     mK = TT*mP * ZZ' * mDinv; %Kalman Gain 
%     
%     amL(:, :, i) = eye(np) - mK * ZZ;
%     
%     % updating
%     alpha = alpha + mK * me(:, i);
%     mP = mP * amL(:, :, i)' + mH2;
%     
% end
% 
% 
% %%--- simulation smoother ---%%
% 
% i = ns;
% while i >= 1
%     mC = mH2 - mH2 * mU * mH2;
%     mC = (mC + mC')/2;
%     [mCc, fl] = chol(mC, 'lower');
%     drC = rcond(mC);
%     if fl > 0
%         mCc = eye(np) * 0.01;
%         mCinv = eye(np) * 10^4;
%     elseif isnan(drC) || (drC < eps*10^2)
%         mCinv = eye(np) * 10^4;
%     else
%         mCinv = inv(mC);
%     end
%     
%     veps = mCc * randn(np, 1);
%     meta(:, i) = mH2 * vr + veps;
%     mV = mH2 * mU * amL(:, :, i);
% 
%     vr = amZ(:,:,i)' * amDinv(:,:,i) * me(:, i) ...
%        + amL(:,:,i)' * vr - mV' * mCinv * veps;
%     mU = amZ(:,:,i)' * amDinv(:,:,i) * amZ(:,:,i) ...
%        + amL(:,:,i)' * mU * amL(:,:,i) + mV' * mCinv * mV;
% 
%     i = i - 1;
% end
% 
% mC = mH02 - mH02 * mU * mH02;
% mC = (mC + mC')/2;
% [mCc, fl] = chol(mC, 'lower');
% if fl > 0
%     mCc = eye(np) * 0.07;
% end
% veta0 = mH02 * vr + mCc *  randn(np, 1);
% 
% malpha = zeros(np, ns);
% malpha(:, 1) = va0 + veta0;
% for i = 1 : ns-1
%     malpha(:, i+1) = malpha(:, i) + meta(:, i);
% end
