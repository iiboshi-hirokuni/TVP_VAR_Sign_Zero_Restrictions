
set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

              svar_type = 'Caldala_Kamp'  % Caldara and Kamp (2017)
%           svar_type = 'Caldala_Kamps_tax'  % Caldara and Kamp (2017)
%             svar_type = 'Mountford_Uhlig'    % Mountford and Uhlig (2009)

ns = m_ns;  
nk = m_nk; 
nl = m_nl; 

us_bc=csvread('./data/US_BC.csv',1);
ti = 1952:0.25:1952+(ns-nl-1)/4;

a0 = load([ char(filename) '/a0_',char(policy_type),'-' ,char(svar_type),'.xls']);

% psi_y = -a0(:,1)./a0(:,4)/ratio_g(1);
% psi_c = -a0(:,2)./a0(:,4)/ratio_g(2);
% psi_i = -a0(:,5)./a0(:,4)/ratio_g(4);
psi_y = -a0(:,1)./a0(:,3);
psi_r = -a0(:,4)./a0(:,3);
psi_g = -a0(:,2)./a0(:,3);
psi_pi = -a0(:,5)./a0(:,3);


figure('File','A0_para_gov_plot');
 h=area(ti, [ -100*ones(ns-nl,1), 500*us_bc(nl+1:end,2)] , 'LineStyle','non') ;
     set(h(1),'FaceColor',[1 1 1])  
     set(h(2),'FaceColor',[0.5,1,1])  % blue

hold on
l1=plot(ti,psi_y(nl+1:end),'k:','LineWidth',2);
l2= plot(ti,psi_pi(nl+1:end),'r-','LineWidth',2);
l3= plot(ti,psi_g(nl+1:end),'k-','LineWidth',1.5);
l4= plot(ti,psi_r(nl+1:end),'b--','LineWidth',2);
hold off
legend([l1,l2,l3,l4],{'\psi_y','\psi_{\pi}','\psi_{\tau}','\psi_r'},'FontSize',12);
ylim([-5 2]);


%%

% figure('File','A0_psi_y_plot');
%  h=area(ti, [ -100*ones(ns-nl,1), 500*us_bc(nl+1:end,2)] , 'LineStyle','non') ;
%      set(h(1),'FaceColor',[1 1 1])  ;
%      set(h(2),'FaceColor',[0.5,1,1]);  % blue
% 
% hold on
% l1=plot(ti,psi_y(nl+1:end),'k:','LineWidth',2);
% 
% hold off
% legend([l1],{'\psi_y'},'FontSize',12)
% ylim([-2 2])
