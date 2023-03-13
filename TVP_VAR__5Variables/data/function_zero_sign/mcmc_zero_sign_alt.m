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

function mcmc_zero_sign(nsim, nburn, jj,m,div_k)     % Aug/10/2013   nburnを追加

% global m.my m_asvar m_nl m_ns m_nk m_fli m_flSb m_nimp m_flfi ...
%        m_iseed m_dvb0 m_dVb0 m_dva0 m_dVa0 m_dvh0 m_dVh0 m_k;
%    
% global policy_type svar_type; 

method_type = m.method_type;
policy_type = m.policy_type;
svar_type   = m.svar_type;
   
tic;

 nthin=1;

%%--- set default options ---%%

if isempty(m.fli) == 1
    m.fli = 0;
end
if isempty(m.flSb) == 1
    m.flSb = 1;
end
if isempty(m.nimp) == 1
    m.nimp = 12 + 1;
end
 if isempty(m.flfi) == 1
     m.flfi = 1;
end
if isempty(m.iseed) == 1
    m.iseed = jj*100;
end

rand('state', m.iseed);



%%--- set variables ---%%

nimp = m.nimp;
ns = m.ns;  % # of time periods
nk = m.nk ; % # of series
nl = m.nl;  % # of lags
nb = nk * (nk*nl + m.fli);  % # of coefficients in beta
na = nk * (nk-1) / 2;       % # of parameters in a

if m.fli == 1
    vym = zeros(1, nk);
else
    vym = zeros(1, nk);
end
m.my = m.my - ones(ns, 1) * vym;

myh = zeros(ns, nk);
mya = zeros(ns, nk);
amX = zeros(nk, nb, ns);
amXh = zeros(nk, na, ns);
amG2 = zeros(nk, nk, ns);
mai = zeros(ns, na);
for i = nl+1 : ns
    amX(:, :, i) = fXt(m.my(i-nl:i-1, :), m.fli);
end

mb = zeros(ns, nb);
ma = zeros(ns, na);
mh = zeros(ns, nk);

mSigb = eye(nb) * 0.01;
mSiga = eye(na) * 0.01;
mSigh = eye(nk) * 0.01;

vidb = 1 : nb;
if m.fli == 1
    vidi = (0 : nk-1) * (nk*nl+1) + 1;
	vidb(vidi) = [];
end
[v1, v2] = find(triu(reshape(1:nk^2, nk, nk)', 1));
vida = (v1-1)*nk + v2;

%%--- prior ---%%

if isempty(m.dvb0) == 1
  if m.flSb == 1
    m.dvb0 = 25;          % Sigma ~ IW(vb0, I*Vb0)
    m.dVb0 = 1e-4;
  else
    m.dvb0 = 40;          % sigb_i^2 ~ IG(va0/2, Va0/2) 
    m.dVb0 = 2*1e-4;
  end
elseif m.flSb == 0
    m.dvb0 = m.dvb0*2;
    m.dVb0 = m.dVb0*2;
end   
if isempty(m.dva0) == 1
  m.dva0 = 8;             % siga_i^2 ~ IG(va0/2, Va0/2)
  m.dVa0 = 2*1e-4;    
end
if isempty(m.dvh0) == 1
  m.dvh0 = 8;             % sigh_i^2 ~ IG(vh0/2, Vh0/2)
  m.dVh0 = 2*1e-4;    
end

vb0 = zeros(nb, 1);       % b_1 ~ N(b0, Sb0)
mSb0 = eye(nb) * 10;
va0 = zeros(na, 1);       % a_1 ~ N(a0, Sa0)
mSa0 = eye(na) * 10;
vh0 = zeros(nk, 1);       % h_1 ~ N(h0, Sh0)
mSh0 = eye(nk) * 50;

mS0 = eye(nb) * m.dVb0;
dnub = m.dvb0 + ns - nl - 1;
dnua = m.dva0 + ns - nl - 1;
dnuh = m.dvh0 + ns - nl - 1;

    
%%--- set sampling option ---%%

% nburn = 0.1 * nsim;         % burn-in period  % Aug/10/2013 停止
npmt = 6;                   % # of parameter to store
msamp    = zeros(nsim, npmt);  % sample box
msamph   = zeros(ns, nk);
msamphs  = zeros(ns, nk);
msampa   = zeros(ns, na);
msampas  = zeros(ns, na);
msampai  = zeros(ns, na);
msampais = zeros(ns, na);

 msampb = zeros(ns, length(vidb));
 msampbs = zeros(ns, length(vidb));

if m.fli == 1
    msampi  = zeros(ns, nk);
    msampis = zeros(ns, nk);
end
if m.flfi == 1
  
else  % % Dec/30/2014 追加
    mimpm_temp = zeros(m.nimp, nk^2,ns);
    mimpms_temp = zeros(m.nimp, nk^2,ns);
   
end

n_IRF = zeros(ns,1);
nK = floor(m.ns/30)-1;      % # of blocks for sampling h
F          = zeros(nk*nl,nk*nl);

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
     = ssmooth(m.my(nl+1:end,:), amX(:,:,nl+1:end), ...
               amG2(:,:,nl+1:end), mSigb, vb0, mSb0)';
    
%       y_temp= m.my(nl+1:end,:)-m.my(nl:end-1,:);  %% 2018/8/21追加
%     
%       mb(nl+1:end, :)  = ssmooth(y_temp, amX(:,:,nl+1:end), ...
%                        amG2(:,:,nl+1:end), mSigb, vb0, mSb0)';  %% 2018/8/21追加
%                    
%           prior_b=reshape([eye(nk), zeros(nk,nk*(nl-1))]',1,nk*nk*nl);
%           
%      mb(nl+1:end, :) = mb(nl+1:end, :) + ones(size(msampb,1)-nl,1).*prior_b;  %% 2018/8/21追加
     

    
  %%--- sampling a ---%%
    
    for i = nl+1 : ns
       myh(i, :) = m.my(i, :) - mb(i, :) * amX(:, :, i)';
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
    
    if m.flSb == 1
      mSb = inv(mS0 + mdif'*mdif);
      mSb = (mSb + mSb')/2;
      [mL, p] = chol(mSb, 'lower');
      if p > 0
        mSb = diag(diag(mSb));
      end
      mSigb = inv(wishrnd(mSb, dnub));
      mSigb = 1*(mSigb + mSigb')/2;
    else
      vSb = m.dVb0 + sum(mdif.^2);
      mSigb = diag(1 ./ gamrnd(dnub/2, 2./vSb));
    end
    
    vSa = m.dVa0 + sum(diff(ma(nl+1:end, :)).^2);
    mSiga = diag(1 ./ gamrnd(dnua/2, 2./vSa));
    
    vSh = m.dVh0 + sum(diff(mh(nl+1:end, :)).^2);
    mSigh = diag(1 ./ gamrnd(dnuh/2, 2./vSh));


%%--- storing sample ---%%
        if mod(m_k, div_k )==0 
             disp([num2str(jj) char('-th CPU :') num2str(m_k), char('-th iteration') ]);
        end        

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
        msampbs = msampbs + mb(:, vidb).^2;
        
        if m.fli == 1
            msampi  = msampi + mb(:, vidi);
            msampis = msampis + mb(:, vidi).^2;
        end
        if m.flfi == 1        
      %%--- impulse response ---%%      
        else
            
%             if mod(m_k, 1 )==0     % Aug/10/2013　追加              
             
        % check of stationarity
           F(1:nk,1:nk*nl)= reshape( mean(mb(nl+1:end,:),1), nk, nk*nl);
               eigen           = max(eig(F));
               largeeig     = abs(eigen);

          if   largeeig < 1   % stationarity
               
               [mimpm_temp2, sign_chk] = ...
                   impulse_zero_sign_v2(nl, m.nimp, mb(:, vidb), ma, mh,...
                                        method_type, policy_type,svar_type);
%               
                for t = 1:ns
                    if sign_chk(t)==1  
                       mimpm_temp(:,:,t) = mimpm_temp(:,:,t) + mimpm_temp2(:,:,t) ;
                       mimpms_temp(:,:,t) = mimpms_temp(:,:,t) + mimpm_temp2(:,:,t).^2 ;
                       n_IRF(t)=n_IRF(t)+1;
                    end
                end                
              
          end                     
            
            
        end
               
    if mod(m_k, div_k) == 0       % print counter
%         fprintf('%i \n', m_k, '-th iterations');
        
     disp(['Accept Rates = ' num2str(100*mean(n_IRF(nl+1:ns)/m_k)), ' (%)'  ] )

        h =figure(999);
         plot(n_IRF(nl+1:ns)/m_k);
%         plot(n_IRF/m_k);
        title('Accept Rates of Sign Restrictons in each period')       
        pause(0.05);
    end
    
    if mod(m_k, div_k*100) == 0   
%          save_impulse;
    end
    
   end  
end


%%--------------- S A M P L I N G   E N D ----------------%%

%%--- output result ---%%

iBm = min([500, nsim/2]);   % bandwidth
iacf = iBm;

aspar = char('sb1  ', 'sb2', 'sa1', 'sa2', 'sh1', 'sh2');
aspar2 = char('  s_{b1}', '  s_{b2}', '  s_{a1}', ...
              '  s_{a2}', '  s_{h1}', '  s_{h2}');
    
est_date = datestr(date);   
result_name = ['./output/result', num2str(nsim), '_', est_date ,policy_type, svar_type, '.txt'];          
fileID = fopen(result_name,'w');
   
fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'------------------------------------');
fprintf(fileID,'\nParameter   Mean      Stdev       ');
fprintf(fileID,'95%%U       95%%L    Geweke     Inef.');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'------------------------------------\n');

if nsim >= 50
msamp = sqrt(msamp);
for i = 1 : npmt
    vsamp = msamp(:, i);
    vsamp_s = sort(vsamp);
fprintf(fileID,'%s %10.4f %10.4f %10.4f %10.4f %9.3f %9.2f\n',...
        aspar(i, :), ...
        [mean(vsamp), std(vsamp), ...
         vsamp_s(floor(nsim*[0.1;0.9]))'], ...
          fGeweke(vsamp, iBm), ...
         ftsvar(vsamp, iBm)/var(vsamp));
end          
end

fprintf(fileID,'-----------------------------------');
fprintf(fileID,'-----------------------------------');
fprintf(fileID,'\nTVP-VAR model (Lag = %i', nl);
fprintf(fileID,')\nIteration: %i', nsim);
if m.flSb == 0
  fprintf(fileID,'\nSigma(b): Diagonal');
end

fclose(fileID);

%%--- output graphs ---%%

%% parameters %%

vacf = zeros(iacf, 1);
figure(101)
for i = 1 : npmt
    for j = 1 : iacf
        macf = corrcoef(msamp(j+1:end, i), ...
                           msamp(1:end-j, i));
        vacf(j) = macf(2, 1);
    end
    subplot(3, npmt, i)        
    sysh = stem(vacf);              % autocorrelation
    set(sysh, 'Marker', 'none')
    axis([0 iacf -1 1])
    title(aspar2(i, :))
    subplot(3, npmt, npmt+i);
    plot(msamp(:, i))               % sample path
    title(aspar2(i, :))
    vax = axis;
    axis([0 nsim vax(3:4)])
    subplot(3, npmt, npmt*2+i)
    hist(msamp(:, i), 15)           % posterior density
    title(aspar2(i, :))
end

%% draw h %%
 
 msamph = msamph / nsim;   % posterior mean
 msamphs = sqrt(msamphs/nsim - msamph.^2);
                          % posterior standard deviation  

if m.fli == 1
    m.my = m.my + ones(ns, 1) * vym;
end

figure(102)
for i = 1 : nk
    subplot(2, nk, i);
    plot(m.my(nl+1:end, i))
    vax = axis;
    axis([0 ns-nl vax(3:4)])
    if vax(3) * vax(4) < 0
        line([0, ns], [0, 0], 'Color', ones(1, 3)*0.6)
    end
    if i == 1
      title(['Data: ', char(m.asvar(i))], ... 
            'interpreter', 'latex')
    else
      title(char(m.asvar(i)), 'interpreter', 'latex')
    end
end
for i = 1 : nk
    subplot(2, nk, i+nk);
    plot(exp(msamph(nl+1:end, i)))
    hold on
    plot(exp(msamph(nl+1:end, i) - msamphs(nl+1:end, i)), 'r:')
    plot(exp(msamph(nl+1:end, i) + msamphs(nl+1:end, i)), 'r:')
    hold off
    vax = axis;
    axis([1 ns-nl vax(3:4)])
    if i == 1
      legend('Posterior mean', '1SD bands')
    else
      title(char(m.asvar(i)), 'interpreter', 'latex')        
    end
end

mout = [msamph msamphs];
save(['./output/tvpvar_vol_' num2str(jj) '.xls'], 'mout', '-ascii', '-tabs');

%% draw b %%
  msampb = msampb / nsim;   % posterior mean
  msampbs = sqrt(msampbs/nsim - msampb.^2);
                         % posterior standard deviation  
  mout = [msampb msampbs];
  save(['./output/tvpvar_b_' num2str(jj) '.xls'], 'mout', '-ascii', '-tabs');

%% draw a %%
 
 msampa = msampa / nsim;   % posterior mean
 msampas = sqrt(msampas/nsim - msampa.^2);
                          % posterior standard deviation  


mout = [msampa msampas];
save(['./output/tvpvar_a_' num2str(jj) '.xls'], 'mout', '-ascii', '-tabs');

%% draw a-inverse %%

 msampai = msampai / nsim;   % posterior mean
 msampais = sqrt(msampais/nsim - msampai.^2);
                          % posterior standard deviation  


mout = [msampa msampas];
save(['./output/tvpvar_ai_' num2str(jj) '.xls'], 'mout', '-ascii', '-tabs');

if m.fli == 1
    
  %% draw intercept %%

  msampi = msampi / nsim;   % posterior mean
  msampis = sqrt(msampis/nsim - msampi.^2);
                          % posterior standard deviation  

  mout = [msampi msampis];
  save( ['./output/tvpvar_int_' num2str(jj) '.xls'], 'mout', '-ascii', '-tabs');
end

%% save impulse response %%

if m.flfi == 1
    mimpm = impulse(nl, m.nimp, msampb, msampa,...
                    msamph);               
           
else
    
    mimpm = [];   
    mimpms = []; 
    for i = 1:size(mimpm_temp,3)
       
        mimpm = [mimpm; mimpm_temp(:,:,i) ] ;
        mimpms = [mimpms;  mimpms_temp(:,:,i) ] ;
        
    end 
    
                      
%     change_shock_size;  %  財政ショックのサイズの基準化   
       nimp = size(mimpm, 1) / ns;
       
     save(['./output/tvpvar_imps_' num2str(jj) '.xls'], 'mimpms',  '-ascii', '-tabs');

end

save(['./output/tvpvar_imp_' num2str(jj) '.xls'], 'mimpm',  '-ascii', '-tabs');

% save_name_para = ['./output/para_save_', num2str(jj) ,'_' , num2str(nsim), '_',...
%                     est_date ,policy_type, svar_type, '.mat'];
save_name_para = ['./output/para_save_', num2str(jj) ,  '.mat'];                
save(save_name_para, 'nl','nimp', 'nsim', 'msampb', 'msampa', 'msamph','n_IRF' );  



fprintf('\n\nRanseed: %i', m.iseed);
fprintf('\nTime: %.2f', toc);
fprintf('\n\n')



