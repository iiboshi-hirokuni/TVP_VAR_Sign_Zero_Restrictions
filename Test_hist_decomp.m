clear all
clc

warning('off','all');

addpath('./data')
addpath('./function')
addpath('./function_zero_sign')
addpath('./output')

load('output/para_save200_20-Feb-2015bench_markboth_zero_sign.mat'); 

my=xlsread('./data/DATA_6var_closed1952');

asvar = {'g';'y';'cons';'dbt';'\pi';'int'};    % variable names
nlag = 4;                   % lags

setvar('data', my, asvar, nlag);  % set data

P = eye(6); 
sign_P = ones(6);

 [decomp_smooth  amOmsq ] = hist_decomp(my, msampa, msampb/nsim, msamph, P, sign_P);
 
 title_data = {'gov' 'y' 'cons' 'debt' 'price' 'ffr'  };
 
 figure(1)
 for i = 1:6 
   subplot(3,2,i)
      plot(my(:,i),'b--')
   hold on
      plot(decomp_smooth(:,i),'r-' )
   hold off
   title( char(title_data(i)) )
   legend('actual','fitted')
 end  
 
  amOmsq
   
