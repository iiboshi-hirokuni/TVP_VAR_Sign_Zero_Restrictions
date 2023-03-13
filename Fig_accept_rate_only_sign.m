
% addpath('../Waseda_results/output_with_zero/')

nsim = 500;
ncores = 4;
n_IRF_total = zeros(ns,1);

for i = 1:ncores  
   
%    load(['../Waseda_results/output_with_zero/para_save_', num2str(i) ,  '.mat']);
   load(['./save_output/para_save_', num2str(i) ,  '.mat']);
   n_IRF_total =  n_IRF_total+ n_IRF;    
   
end   


 fig1= figure(1000);
 ti=1953:0.25:1953+(ns-nl-1)/4;
 plot(ti, n_IRF_total(nl+1:ns)/nsim/ncores*100);
%   title('Accepted Rates for Zero and Sign Restictions')
  title('Accepted Rates for only Sign Restictions')
  ylabel("%")
  ylim([0 25])

  est_date = datestr(date);   
   name = ['./Fig/Accepted Rates_',est_date];
    saveas(fig1,name,'fig');    
  