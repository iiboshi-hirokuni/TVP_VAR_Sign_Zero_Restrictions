%%--------------------------------------------------------%%
%%                     tvpvar_ex1.m                       %%
%%--------------------------------------------------------%%

%%
%%  MCMC estimation for Time-Varying Parameter VAR model
%%  with stochastic volatility
%%
%%  tvpvar_ex*.m illustrates MCMC estimation
%%  using TVP-VAR Package
%%  (Data: "tvpvar_ex.xls/.mat")
%%

clear all;
close all;
clc

warning('off','all')


global method_type policy_type svar_type;

addpath('./data')
addpath('./function')
addpath('./function_zero_sign')
addpath('./output')

% my=xlsread('./data/data_7vars_qdt','sheet1','B2:H266');
my=csvread('./data/data_5vars_ldt.csv',1,1);

%% =====================================================
   method_type = 'Arias_et_al'

   policy_type = 'bench_mark' 
%  policy_type = 'Choleski' 

         svar_type = 'Caldala_Kamp'  % Caldara and Kamp (2017)
%             svar_type = 'Caldala_Kamps_tax' 
%           svar_type = 'Mountford_Uhlig'    % Mountford and Uhlig (2009)

%% =====================================================

asvar = {'y';'tax';'g';'int';'\pi'};    % variable names
nlag = 2;                   % lags

setvar('data', my, asvar, nlag);  % set data

%setvar('ranseed', 5);       % set ranseed
setvar('intercept', 0);     % set time-varying intercept
setvar('SigB', 1);          % set digaonal for Sig_beta
setvar('impulse', 20);      % set maximum length of impulse
                            % (smaller, calculation faster)
                            
setvar('prior','b',20, 1e-4);   % set prior of var of b_t 
   setvar('prior','a',10, 1*1e-4);  % set prior of var of a_t 
  setvar('prior','h',50, 1e-4);  % set prior of var of h_t 

setvar('fastimp', 0);       % fast computing of response
                         % 0:各サンプルからインパルスを計算し、インパルスの平均値を計算
                         % 1: サンプルの平均値からインパルスを計算　%  Aug/10/2013 追加                         
 change_global;

% plot_data;   
pause(0.5);

%% MCMC
% tic
% 
    ncores = 4;
%    delete(gcp('nocreate')) %   
%    parpool('local', ncores)
%   
 jj = 1;
%  parfor jj = 1:ncores
         mcmc_zero_sign(100, 10, jj,m, 10);   %　thining(間引き数)はなし  
%   end                 
%  
% toc
%%

integrate_save_files       

%%      
% %%   Check of Convergence
% 
% F          = zeros(m_nk*m_nl,m_nk*m_nl);
% load('./output/para_save_1.mat');
% 
% disp(' Coeffi_s of reduced VAR = ')
% coef_B = reshape( mean(msampb(nl+2:end,:),1), m_nk*m_nl, m_nk)'
% F(1:m_nk,1:m_nk*m_nl)= coef_B;
% 
%      eigen           = max(eig(F(:,1:m_nk*m_nl)));
%      eigen           = max(eigen);
%      largeeig     = abs(eigen);
%      disp([' eig = ', num2str(largeeig)    ]);
% 
%      
%%  Make Graph 
filename=('./output/'); 
%    drawimp_3d;  
%   pause(5)
% %  
       drawimp_68band;
% %     pause(5)
%     
%      draw_Fiscal_Mult
%
if svar_type == 'Caldala_Kamp'
    draw_a0
elseif svar_type == 'Caldala_Kamps_tax'
    draw_a0_tax
end

%%  save result 
%  filename=('./output/'); 
% mimpm = load([ filename '/tvpvar_imp_',char(policy_type),'-' ,char(svar_type),'.xls']);     % posterior means of IRF
% mimpms = load( [ filename '/tvpvar_imps_',char(policy_type),'-' ,char(svar_type),'.xls' ] );   % Standard Deviation of IRF
%  accept_rate = load( [ filename '/accept_rate_',char(policy_type),'-' ,char(svar_type),'.xls' ] );
%  save(['./save_output/tvpvar_imp_',char(policy_type),'-' ,char(svar_type),'.xls'],...
%                         'mimpm',  '-ascii', '-tabs');
% save(['./save_output/tvpvar_imps_',char(policy_type),'-' ,char(svar_type),'.xls'],...
%                         'mimpms',  '-ascii', '-tabs');
%  save(['./save_output/accept_rate_',char(policy_type),'-' ,char(svar_type),'.xls'],...
%                         'accept_rate',  '-ascii', '-tabs');  
%  save(['./save_output/a0_',char(policy_type),'-' ,char(svar_type),'.xls'],...
%                         'a0_mean',  '-ascii', '-tabs');                      

            
                    
