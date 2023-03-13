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

global method_type policy_type svar_type;

addpath('./data')
addpath('./function')
addpath('./function_zero_sign')
addpath('./output')

% my=xlsread('./data/data_7vars_qdt','sheet1','B2:H266');
my=csvread('./data/data_5vars_ldt.csv',1,1);


%% =====================================================

   method_type = 'Arias_et_al'
%  method_type= 'only_sign'

policy_type = 'bench_mark' 
%      policy_type = 'unanticipated_fiscal_policy'
%    policy_type = 'anticipated_fiscal_policy'

%          svar_type = 'Caldala_Kamp'  % Caldara and Kamp (2017)
         svar_type = 'Caldala_Kamps_tax' 
%      svar_type = 'Mountford_Uhlig'    % Mountford and Uhlig (2009)

%% =====================================================

asvar = {'y';'tax';'g';'int';'\pi'};    % variable names
nlag = 1;                   % lags

setvar('data', my, asvar, nlag);  % set data

%setvar('ranseed', 5);       % set ranseed
setvar('intercept', 1);     % set time-varying intercept
setvar('SigB', 0);          % set digaonal for Sig_beta
setvar('impulse', 80);      % set maximum length of impulse
                            % (smaller, calculation faster)

setvar('fastimp', 0);       % fast computing of response
                         % 0:各サンプルからインパルスを計算し、インパルスの平均値を計算
                         % 1: サンプルの平均値からインパルスを計算　%  Aug/10/2013 追加                         
change_global;
  plot_data;    
  
%%
% filename=('./save_output_80/'); 
 filename=('./save_output/'); 
 
  drawimp_68band; 

  switch svar_type
      case  'Caldala_Kamp'
             draw_a0         
%               draw_Fiscal_Mult; 
      case  'Caldala_Kamps_tax'
             draw_a0_tax
%             draw_Fiscal_Mult_tax;  
  end      
        pause(0.5)  
        
        
%          drawimp_3d;  
%   pause(0.5)      
    

%%  
  filename1=('./save_output/'); 
   filename2=('./save_output/');
% % filename2=('./output/');
% %  
    Fig_accept_rate;
   draw_Fiscal_Mult_compare_model;



