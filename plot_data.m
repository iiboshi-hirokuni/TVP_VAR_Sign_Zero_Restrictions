% ==============================================
% 
%     Plot data
% 
% ==============================================

global m_ns m_nk m_nl m_asvar m_vpd policy_type;

type = 2;  % 2: Recessions 
           % 3: Active Monetary Policy regime
           % 4: Active Fiscal Policy 
           % 5: Passive Fiscal Policy 
           % 6: Political Regime 
           
% load('para_save.mat');

us_bc=csvread('./data/US_BC.csv',1);

figure(1)
y = 1952:0.25:1952+(m_ns-1)/4;
 ti = 1952+(m_nl+4+1)/4:0.25:1952+(m_ns)/4;
 
for i = 1:m_nk
  subplot(3,2,i)
  
  h=area(ti, [ -100*ones(m_ns-(m_nl+4),1) 500*us_bc(m_nl+4+1:end,type)] , 'LineStyle','non') ;
      set(h(1),'FaceColor',[1 1 1])  
    if type == 2  
      set(h(2),'FaceColor',[1,0.5,1])  % blue
    elseif type == 3  
       set(h(2),'FaceColor',[0.5,1,0.5]) % Green
    elseif type == 5  
        set(h(2),'FaceColor',[1, 0.75, 0.5])  % Yellow            
    elseif type == 6  
         set(h(2),'FaceColor',[1,0.5,1])  % blue
    end
    
    l2 = h(2);
    
   hold on
      l1=plot(y(m_nl+4+1:end),my(m_nl+4+1:end,i),'LineStyle','-','Color','b',...
        'LineWidth',2.5);
   hold off  
    
  title([ char(m_asvar(i)) ],'FontSize',13,'FontWeight','bold');  
     xlim([1950 2020])
     if i >= 6
          ylim([-0.0 0.2]) 
     else    
         ylim([-0.2 0.2])
     end    
  if i == 6
      if type==6
        legend([l1,l2], {'data','Repubilic'} )
      elseif type ==2
        legend([l1,l2], {'data','Recessions'},'FontSize',10 )
      end    
 end    
end