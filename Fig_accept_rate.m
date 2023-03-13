
% addpath('../Waseda_results/output_with_zero/')
ns = m.ns;
nl = m.nl;
% nsim = 100;
% ncores = 16;


for j = 1:2
   
  if j == 1  
       filename =  filename1;
          svar_type = 'Caldala_Kamp'
  elseif  j== 2
      filename =  filename2;
%       svar_type = 'Mountford_Uhlig'
        svar_type = 'Caldala_Kamps_tax'
  elseif  j== 3
      filename =  filename3;  
  elseif  j== 4
      filename =  filename4;
  end    
  
   
  accept(:,j)=  load([ filename '/accept_rate_',char(policy_type),'-' ,char(svar_type),'.xls']);

end

 fig1= figure(1000);
 ti=(1952+nl/4):0.25:(1952+(ns-2)/4); 
     L1=plot(ti,accept(:,1),'LineStyle','-','Color','r', 'LineWidth',1.5); 
 hold on
     L2=plot(ti, accept(:,2),'LineStyle','--','Color','b', 'LineWidth',1.5); 
 hold off
%   title('Accepted Rates for Zero and Sign Restictions')
  title({'Accepted Rates of Zeros and Sign Restrictions'},'FontSize',12)
%   legend('w/ Zero, All shocks','w/ Zero, FP+BC shocks', 'w/ Zero, FP shocks','w/o Zero, All shocks' )
%   legend('Model V2','Model V2 prime','Model V3','Model V3 prime')
   legend([L1,L2],{'Gov. Spending  Rule','Tax Cut Rule'},'FontSize',11)
  ylabel({"[ % ]"}, 'FontSize',11)
  xlim([1950 2020])
   ylim([0 25])

  est_date = datestr(date);   
   name = ['./Fig/Accepted Rates_',est_date];
    saveas(fig1,name,'fig');    
  