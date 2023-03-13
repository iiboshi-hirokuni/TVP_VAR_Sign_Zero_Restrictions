%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  [] = mcmc(nsim)
%%
%%  "mcmc" implements MCMC estimation for TVP-VAR model
%%
%%  [input]
%%      nsim:  # of MCMC iterations
%%

function mcmc(nsim, nburn, nthin)     % Aug/10/2013   nburnÇí«â¡

global m_my m_asvar m_nl m_ns m_nk m_fli m_flSb m_nimp m_flfi ...
       m_iseed m_dvb0 m_dVb0 m_dva0 m_dVa0 m_dvh0 m_dVh0 m_k
   
tic;

%%--- set default options ---%%

if isempty(m_fli) == 1
    m_fli = 0;
end
if isempty(m_flSb) == 1
    m_flSb = 1;
end
if isempty(m_nimp) == 1
    m_nimp = 12 + 1;
end
 if isempty(m_flfi) == 1
     m_flfi = 1;
end
if isempty(m_iseed) == 1
    m_iseed = 1;
end

rand('state', m_iseed);



%%--- set variables ---%%

ns = m_ns;  % # of time periods
nk = m_nk;  % # of series
nl = m_nl;  % # of lags
nb = nk * (nk*nl + m_fli);  % # of coefficients in beta
na = nk * (nk-1) / 2;       % # of parameters in a

if m_fli == 1
    vym = zeros(1, nk);
else
    vym = mean(m_my);
end
m_my = m_my - ones(ns, 1) * vym;

myh = zeros(ns, nk);
mya = zeros(ns, nk);
amX = zeros(nk, nb, ns);
amXh = zeros(nk, na, ns);
amG2 = zeros(nk, nk, ns);
mai = zeros(ns, na);
for i = nl+1 : ns
    amX(:, :, i) = fXt(m_my(i-nl:i-1, :), m_fli);
end

mb = zeros(ns, nb);
ma = zeros(ns, na);
mh = zeros(ns, nk);

mSigb = eye(nb) * 0.01;
mSiga = eye(na) * 0.01;
mSigh = eye(nk) * 0.01;

vidb = 1 : nb;
if m_fli == 1
    vidi = (0 : nk-1) * (nk*nl+1) + 1;
	vidb(vidi) = [];
end
[v1, v2] = find(triu(reshape(1:nk^2, nk, nk)', 1));
vida = (v1-1)*nk + v2;

%%--- prior ---%%

if isempty(m_dvb0) == 1
  if m_flSb == 1
    m_dvb0 = 25;          % Sigma ~ IW(vb0, I*Vb0)
    m_dVb0 = 1e-4;
  else
    m_dvb0 = 40;          % sigb_i^2 ~ IG(va0/2, Va0/2) 
    m_dVb0 = 2*1e-4;
  end
elseif m_flSb == 0
    m_dvb0 = m_dvb0*2;
    m_dVb0 = m_dVb0*2;
end   
if isempty(m_dva0) == 1
  m_dva0 = 8;             % siga_i^2 ~ IG(va0/2, Va0/2)
  m_dVa0 = 2*1e-4;    
end
if isempty(m_dvh0) == 1
  m_dvh0 = 8;             % sigh_i^2 ~ IG(vh0/2, Vh0/2)
  m_dVh0 = 2*1e-4;    
end

vb0 = zeros(nb, 1);       % b_1 ~ N(b0, Sb0)
mSb0 = eye(nb) * 10;
va0 = zeros(na, 1);       % a_1 ~ N(a0, Sa0)
mSa0 = eye(na) * 10;
vh0 = zeros(nk, 1);       % h_1 ~ N(h0, Sh0)
mSh0 = eye(nk) * 50;

mS0 = eye(nb) * m_dVb0;
dnub = m_dvb0 + ns - nl - 1;
dnua = m_dva0 + ns - nl - 1;
dnuh = m_dvh0 + ns - nl - 1;

    
%%--- set sampling option ---%%

% nburn = 0.1 * nsim;         % burn-in period  % Aug/10/2013 í‚é~
npmt = 6;                   % # of parameter to store
msamp    = zeros(nsim, npmt);  % sample box
msamph   = zeros(ns, nk);
msamphs  = zeros(ns, nk);
msampa   = zeros(ns, na);
msampas  = zeros(ns, na);
msampai  = zeros(ns, na);
msampais = zeros(ns, na);

 msampb = zeros(ns, length(vidb));

if m_fli == 1
    msampi  = zeros(ns, nk);
    msampis = zeros(ns, nk);
end
if m_flfi == 1
   
else
    mimpm = zeros(ns*m_nimp, nk^2);
end
nK = floor(m_ns/30)-1;      % # of blocks for sampling h


%%--- MCMC sampling ---%%

fprintf('\nIteration:\n');

%%------------- S A M P L I N G   S T A R T --------------%%

for m_k = -nburn : nsim

  %%--- sampling beta ---%%

    for i = nl+1 : ns
        mAinv = finvm(fAt(ma(i, :), nk));
        amG2(:, :, i) = mAinv * diag(exp(mh(i,:))) * mAinv';
        mai(i, :) = mAinv(vida)';
    end
  
    mb(nl+1:end, :) ...
     = ssmooth(m_my(nl+1:end,:), amX(:,:,nl+1:end), ...
               amG2(:,:,nl+1:end), mSigb, vb0, mSb0)';
    
    
  %%--- sampling a ---%%
    
    for i = nl+1 : ns
       myh(i, :) = m_my(i, :) - mb(i, :) * amX(:, :, i)';
       amXh(:, :, i) = fXh(myh(i, :), nk, na);
       amG2(:, :, i) = diag(exp(mh(i, :)));
    end
  
    ma(nl+1:end, :) ...
     = ssmooth(myh(nl+1:end,:), amXh(:,:,nl+1:end), ...
               amG2(:,:,nl+1:end), mSiga, va0, mSa0)';
  
  %%--- sampling h ---%%

    for i = nl+1 : ns
        mya(i, :) = myh(i, :) * fAt(ma(i, :), nk)';
    end
           
    for i = 1 : nk
        mh(nl+1:end, i) ...
         = svsamp(mya(nl+1:end,i), mh(nl+1:end,i), ...
                  mSigh(i,i), vh0(i), mSh0(i,i), nK);
    end


  %%--- sampling Sigma ---%%
  
    mdif = diff(mb(nl+1:end, :));
    if m_flSb == 1
      mSb = inv(mS0 + mdif'*mdif);
      mSb = (mSb + mSb')/2;
      [mL, p] = chol(mSb, 'lower');
      if p > 0
        mSb = diag(diag(mSb));
      end
      mSigb = inv(wishrnd(mSb, dnub));
      mSigb = (mSigb + mSigb')/2;
    else
      vSb = m_dVb0 + sum(mdif.^2);
      mSigb = diag(1 ./ gamrnd(dnub/2, 2./vSb));
    end
    
    vSa = m_dVa0 + sum(diff(ma(nl+1:end, :)).^2);
    mSiga = diag(1 ./ gamrnd(dnua/2, 2./vSa));
    
    vSh = m_dVh0 + sum(diff(mh(nl+1:end, :)).^2);
    mSigh = diag(1 ./ gamrnd(dnuh/2, 2./vSh));


%%--- storing sample ---%%

    if m_k > 0
        msamp(m_k, :) = [mSigb(1, 1) mSigb(2, 2) ...
                         mSiga(1, 1) mSiga(2, 2) ...
                         mSigh(1, 1) mSigh(2, 2)];

        msamph   = msamph  + mh;
        msamphs  = msamphs + mh.^2;
        msampa   = msampa  + ma;
        msampas  = msampas + ma.^2;
        msampai  = msampai  + mai;
        msampais = msampais + mai.^2;

        msampb = msampb + mb(:, vidb);
        
        if m_fli == 1
            msampi  = msampi + mb(:, vidi);
            msampis = msampis + mb(:, vidi).^2;
        end
        if m_flfi == 1        
      %%--- impulse response ---%%      
        else
            
            if mod(m_k, nthin )==0     % Aug/10/2013Å@í«â¡
                
               fprintf('%i \n', m_k);
               
             
                mimpm = mimpm ...
                  + impulse(nl, m_nimp, mb(:, vidb), ma, mh);
                
              
            end                       % Aug/10/2013Å@í«â¡
            
            
        end
        
    end
        
    if mod(m_k, 100) == 0       % print counter
        fprintf('%i \n', m_k);
    end

end

%%--------------- S A M P L I N G   E N D ----------------%%

%%--- output result ---%%

iBm = min([500, nsim/2]);   % bandwidth
iacf = iBm;

aspar = char('sb1  ', 'sb2', 'sa1', 'sa2', 'sh1', 'sh2');
aspar2 = char('  s_{b1}', '  s_{b2}', '  s_{a1}', ...
              '  s_{a2}', '  s_{h1}', '  s_{h2}');
    
    
% fprintf('\n\n                        [ESTIMATION RESULT]')
% fprintf('\n----------------------------------')
% fprintf('------------------------------------')
% fprintf('\nParameter   Mean      Stdev       ')
% fprintf('95%%U       95%%L    Geweke     Inef.')
% fprintf('\n----------------------------------')
% fprintf('------------------------------------\n')

msamp = sqrt(msamp);
% for i = 1 : npmt
%     vsamp = msamp(:, i);
%     vsamp_s = sort(vsamp);
% fprintf('%s %10.4f %10.4f %10.4f %10.4f %9.3f %9.2f\n',...
%         aspar(i, :), ...
%         [mean(vsamp), std(vsamp), ...
%          vsamp_s(floor(nsim*[0.1;0.9]))'], ...
% %          fGeweke(vsamp, iBm), ...
%          ftsvar(vsamp, iBm)/var(vsamp))
% end          

% fprintf('-----------------------------------')
% fprintf('-----------------------------------')
% fprintf('\nTVP-VAR model (Lag = %i', nl)
% fprintf(')\nIteration: %i', nsim)
% if m_flSb == 0
%   fprintf('\nSigma(b): Diagonal')
% end


%%--- output graphs ---%%

%% parameters %%
% 
% vacf = zeros(iacf, 1);
% for i = 1 : npmt
%     for j = 1 : iacf
%         macf = corrcoef(msamp(j+1:end, i), ...
%                            msamp(1:end-j, i));
%         vacf(j) = macf(2, 1);
%     end
%     subplot(3, npmt, i)        
%     sysh = stem(vacf);              % autocorrelation
%     set(sysh, 'Marker', 'none')
%     axis([0 iacf -1 1])
%     title(aspar2(i, :))
%     subplot(3, npmt, npmt+i);
%     plot(msamp(:, i))               % sample path
%     title(aspar2(i, :))
%     vax = axis;
%     axis([0 nsim vax(3:4)])
%     subplot(3, npmt, npmt*2+i)
%     hist(msamp(:, i), 15)           % posterior density
%     title(aspar2(i, :))
% end
% 
%% draw h %%
 
 msamph = msamph / nsim;   % posterior mean
 msamphs = sqrt(msamphs/nsim - msamph.^2);
                          % posterior standard deviation  

if m_fli == 1
    m_my = m_my + ones(ns, 1) * vym;
end

% figure
% for i = 1 : nk
%     subplot(2, nk, i);
%     plot(m_my(nl+1:end, i))
%     vax = axis;
%     axis([0 ns-nl vax(3:4)])
%     if vax(3) * vax(4) < 0
%         line([0, ns], [0, 0], 'Color', ones(1, 3)*0.6)
%     end
%     if i == 1
%       title(['Data: ', char(m_asvar(i))], ... 
%             'interpreter', 'latex')
%     else
%       title(char(m_asvar(i)), 'interpreter', 'latex')
%     end
% end
% for i = 1 : nk
%     subplot(2, nk, i+nk);
%     plot(exp(msamph(nl+1:end, i)))
%     hold on
%     plot(exp(msamph(nl+1:end, i) - msamphs(nl+1:end, i)), 'r:')
%     plot(exp(msamph(nl+1:end, i) + msamphs(nl+1:end, i)), 'r:')
%     hold off
%     vax = axis;
%     axis([1 ns-nl vax(3:4)])
%     if i == 1
%       title(['SV $\sigma_t^2=\exp(h_t)$: ', ...
%              char(m_asvar(i))], 'interpreter', 'latex')
%       legend('Posterior mean', '1SD bands')
%     else
%       title(char(m_asvar(i)), 'interpreter', 'latex')        
%     end
% end
% 
% mout = [msamph msamphs];
% save('tvpvar_vol.xls', 'mout', '-ascii', '-tabs');
% 
%% draw a %%
 
 msampa = msampa / nsim;   % posterior mean
 msampas = sqrt(msampas/nsim - msampa.^2);
%                           % posterior standard deviation  
% figure
% for i = 1 : na
%     subplot(ceil(na/2), 2, i);
%     plot(msampa(nl+1:end, i))
%     hold on
%     plot(msampa(nl+1:end, i) - msampas(nl+1:end, i), 'r:')
%     plot(msampa(nl+1:end, i) + msampas(nl+1:end, i), 'r:')
%     hold off
%     vax = axis;
%     axis([0 ns-nl vax(3:4)])
%     if vax(3) * vax(4) < 0
%         line([nl-1, ns+1], [0, 0], 'Color', ones(1, 3)*0.6)
%     end
%     title(['$a_{', num2str(i), 't}$'], ...
%           'interpreter', 'latex')
%     if i == 1
%       legend('Posterior mean', '1SD bands')
%     end
% end
% 
% mout = [msampa msampas];
% save('tvpvar_a.xls', 'mout', '-ascii', '-tabs');
% 
% %% draw a-inverse %%
% 
 msampai = msampai / nsim;   % posterior mean
 msampais = sqrt(msampais/nsim - msampai.^2);
%                           % posterior standard deviation  
% figure
% for i = 1 : na
%     subplot(ceil(na/2), 2, i);
%     plot(msampai(nl+1:end, i))
%     hold on
%     plot(msampai(nl+1:end, i) - msampais(nl+1:end, i), 'r:')
%     plot(msampai(nl+1:end, i) + msampais(nl+1:end, i), 'r:')
%     hold off
%     vax = axis;
%     axis([0 ns-nl vax(3:4)])
%     if vax(3) * vax(4) < 0
%         line([nl-1, ns+1], [0, 0], 'Color', ones(1, 3)*0.6)
%     end
%     if i == 1
%       title(['$\tilde{a}_{', num2str(i), ...
%              't}$ ($A_t^{-1}$: ', ...
%              char(m_asvar(fix((vida(i)-1)/nk)+1)), ...
%              '$\to$', ...
%              char(m_asvar(mod(vida(i)-1,nk)+1)), ')'], ...
%              'interpreter', 'latex')
%       legend('Posterior mean', '1SD bands')
%     else
%       title(['$\tilde{a}_{', num2str(i), 't}$ (', ...
%              char(m_asvar(fix((vida(i)-1)/nk)+1)), ...
%              '$\to$', ...
%              char(m_asvar(mod(vida(i)-1,nk)+1)), ')'], ...
%              'interpreter', 'latex')
%     end
% end
% 
% mout = [msampa msampas];
% save('tvpvar_ai.xls', 'mout', '-ascii', '-tabs');
% 
% if m_fli == 1
%     
%   %% draw intercept %%
% 
%   msampi = msampi / nsim;   % posterior mean
%   msampis = sqrt(msampis/nsim - msampi.^2);
%                           % posterior standard deviation  
% 
%   figure
%   for i = 1 : nk
%     subplot(ceil(nk/2), 2, i);
%     plot(msampi(nl+1:end, i))
%     hold on
%     plot(msampi(nl+1:end, i) - msampis(nl+1:end, i), 'r:')
%     plot(msampi(nl+1:end, i) + msampis(nl+1:end, i), 'r:')
%     hold off
%     vax = axis;
%     axis([0 ns-nl vax(3:4)])
%     if vax(3) * vax(4) < 0
%         line([nl-1, ns+1], [0, 0], 'Color', ones(1, 3)*0.6)
%     end
%     if i == 1
%       title(['Intercept ($c_t$): ', char(m_asvar(i))], ...
%              'interpreter', 'latex')
%       legend('Posterior mean', '1SD bands')
%     else
%       title(char(m_asvar(i)), 'interpreter', 'latex')        
%     end
%   end
% 
%   mout = [msampi msampis];
%   save('tvpvar_int.xls', 'mout', '-ascii', '-tabs');
% end
% 
% %% save impulse response %%
% 
if m_flfi == 1
    mimpm = impulse(nl, m_nimp, msampb/nsim, msampa,...
                    msamph);
                
           
else
    %     mimpm = mimpm / nsim;
    mimpm = mimpm / (nsim / nthin) ; % Aug/10/2013Å@í«â¡
end

save('tvpvar_imp.xls', 'mimpm',  '-ascii', '-tabs');

save('para_save.mat', 'nl','m_nimp', 'nsim', 'msampb', 'msampa', 'msamph' );   


fprintf('\n\nRanseed: %i', m_iseed);
fprintf('\nTime: %.2f', toc);
fprintf('\n\n')

