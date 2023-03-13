% =================================================================
% 
%    財政ショックのサイズの基準化 
% 
% =================================================================

% style = 'run_alone';　% 読み込むエクセルファイルを設定して、単独で起動する。

 style = 'sub-routine'; %% mcmc_sign.m の中で　処理を実行

% style = 'save_Yes';  % run_alone'; の処理　をせずに　PCに保存されている mimpデータを修正して
                     % 計算したインパルスをtvpvar_imp.xlsとして保存 

% global  m_nl m_ns m_nk policy_type ;

switch style
  case 'run_alone'    

   file_name = '1000_26-Jan-2015_unanticipated_fiscal_policy'
   file_name_para = '1000_26-Jan-2015unanticipated_fiscal_policy'

   mimpm = load(['./output/tvpvar_imp_', char(file_name) ,'.xls']);     % posterior means of IRF
   mimpms = load(['./output/tvpvar_imps_',char(file_name) ,'.xls']);   % Standard Deviation of IRF
  load(['./output/para_save',char(file_name_para),'.mat']);

%   ns = m_ns;  % # of time periods
%   nk = m_nk;  % # of series
%   nl = m_nl;  % # of lags
%   nimp = m_nimp;
  
case 'sub-routine'   
        nimp = size(mimpm, 1) / ns;
%         ns = ns; 
case 'save_Yes'   
        nimp = size(mimpm, 1) / ns;
%         ns = ns;      
end

switch policy_type
    case 'bench_mark'
        FP_shock = 3;
end        

for i = 1 : nk  % shock, 1:Gov, 2:BC, 6:Monetary Policy 

    for j = 1 : nk    % number of variables of time series 
     
     id = (i-1)*nk + j;
     mimp = reshape(mimpm(:, id), nimp, ns)';   % posterior mean     
     mimps = reshape(mimpms(:, id), nimp, ns)'; % standard dev
     
    %% 財政ショックのサイズの基準化  ===========
   
     if (i == FP_shock)&(j==1) % i = shock(FP), j=variable(Gov)
         FP = 1./ mimp(:,1);
  
         mimp = (FP*ones(1,nimp)).*mimp;  % posterior mean       
         mimps = (FP*ones(1,nimp)).*mimps; % standard dev
     
     elseif (i == FP_shock)&(j>1)
       if j==2
           ratio = 4.77476  ; %  ratio = Y/G 
       elseif j==3
           ratio = 3.017358 ; %  ratio = C/G 
       else
           ratio = 1;
       end           
         mimp = ratio.*(FP*ones(1,nimp)).*mimp;  % posterior mean
         mimps = ratio.*(FP*ones(1,nimp)).*mimps; % standard dev
     
     end
     
    %% =====================================
     
    mimpm(:,id) = reshape(mimp',nimp*ns,1);
    mimpms(:,id) = reshape(mimps',nimp*ns,1);
    
    end
  
    
end 

% switch style
%     case 'run_alone'  
%       save('./output/tvpvar_imp.xls', 'mimpm',  '-ascii', '-tabs');
%       save('./output/tvpvar_imps.xls', 'mimpms',  '-ascii', '-tabs');
%     case 'save_Yes'
%       save('./output/tvpvar_imp.xls', 'mimpm',  '-ascii', '-tabs');
%       save('./output/tvpvar_imps.xls', 'mimpms',  '-ascii', '-tabs'); 
% end
         