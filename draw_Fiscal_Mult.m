%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%

% discount factor
avri = 0.04314;        %average of nominal interest rate
beta = 1/(1+avri); 

        %   Y/G         Tax/G  G/G  int/G  inf/g
ratio_g =[ 4.83099,  1/1,   1/1,  1,     1    ];

%%  インパルス応答の作成

 j = 1;   % 出力したい変数; 1:Gov, 2:Y, 3:cons, 4:debt, 5:pi, 6:int 

global m_ns m_nk m_nl m_asvar m_vpd policy_type;

ns = m_ns;  
nk = m_nk; 
nl = m_nl; 

type = 2;  % 2: Recessions 
           % 3: Active Monetary Policy regime
           % 4: Active Fiscal Policy 
           % 5: Passive Fiscal Policy 
           
% load('para_save.mat');

us_bc=csvread('./data/US_BC.csv',1);

nimp = size(mimpm, 1) / m_ns;
mline = [1 0 0;  0 0 1; 0 .5 0;  0 .7 .7];
vline = {':.', '-', '--', ':'};
% nline = size(vt, 2);

var_seq = { 'Y';'Tax'; 'Gov'; 'R'; 'Pi' };
 
switch  policy_type  
 case 'bench_mark'
    
     graph_start= 1;
     graph_end  = 4;
     idx = 2;    
     switch svar_type
        case 'Mountford_Uhlig'  
             var_Gov = 4;
            shock_gov = 3;
            shock_seq =[1 2 3 4];  % 1:BC shock, 2:MP shock, 3:FP shock, 4:TAX 
            asshock = {'BC';'MP';'Gov';'TAX';' ';'real';'Demand';'non'};
        case 'Caldala_Kamp' 
         var_Gov = 3;
            shock_gov = 1;
            shock_seq =[1 2 3 4];  % 1:FP shock 
            asshock = {'gov';'TAX';' ';'real';'Demand';'non'};
     end    
   
end

time_title={'1960Q1' '1970Q1' '1980Q1' '1990Q1' '2000Q1' '2010Q1' }'; 
time_seq  =[33 73 113 153 193 233]; 

ti = 1952:0.25:1952+(ns-nl-1)/4;
horizon =[1 4  12  20];
hz = 0:1:(nimp-1); 

 
 %  期間による累積の財政乗数の計算
 
 vimp=[3  5  13 21];  % インパルスを合計する期間

cha_imp = {'0.5 years';'1 years';...
           '3 years'; '5 years' };
cha_imp_1 = {'2';'4';'8';'12'; '20' };

save_multp=zeros(size(vimp,2),nk,ns); 
save_multps=zeros(size(vimp,2),nk,ns);
multp=zeros(ns,1);
multps=zeros(ns,1);

for i = 1: size(vimp,2) 
% h_(4000+10*i) = figure(4000+10*i);

   for j = 1 : nk 
     idg = (shock_seq(shock_gov)-1)*nk +  var_Gov;  
    id = (shock_seq(shock_gov)-1)*nk + j;
     mimp_g = reshape(mimpm(:, idg), nimp, ns)';     
    mimp = reshape(mimpm(:, id), nimp, ns)';
    
     mimps_g = reshape(mimpms(:, idg), nimp, ns)';     
    mimps = reshape(mimpms(:, id), nimp, ns)';
    
   
    for k = nl+1:ns
       
       l= vimp(i); %nimp;  % インパルスを合計する期間
%         multp(k)=sum(mimp(k,1:l),2)/sum(mimp_g(k,1:l),2);       
       
        %% 乗数の現在割引価値 
         multp(k)=(sum_pv(mimp(k,1:l),beta)/sum_pv(mimp_g(k,1:l),beta))*ratio_g(j);  
%         multp(k)=mean(mimp(k,1:l),2); 
       
       %  multipliers for variable i = Sum(impluse of variable i) /  Sum(impluse of gov cons + inv)
       
         multps(k)=sum(mimps(k,1:l),2)/sum(mimp_g(k,1:l),2);  
        
%        multps(k)=mean(mimps(k,1:l),2) ;  
      
    end
    
     save_multp(i,j,:) = multp;
     save_multps(i,j,:) = multps;
    
%     for l = nl+1+3:ns-2
%       multp_ma(l)= (0.5*(multp(l-2)+multp(l+2))+sum(multp(l-1:l+1)))/4;  %% 4期間の移動平均
%       
%         multps_ma(l)= (0.5*(multps(l-2)+multps(l+2))+sum(multps(l-1:l+1)))/4;
%     end  

   end
   
end

for j=1:nk
     h_(5000+10*j) = figure(4000+10*j);
      set( h_(5000+10*j),'Position',[20,20,900,600]);
      
  for i = 1:size(vimp,2) 
      
      multp= squeeze(save_multp(:,j,:));
      
     subplot(2,2, i);
%       subplot(nk/3,3, j);
     ti = 1952+(m_nl+4+1)/4:0.25:1952+(m_ns)/4;
%      ti = 1953+(nl+4)/4:0.25:1953+(ns-2)/4;
     
      h=area(ti, [ -100*ones(m_ns-(m_nl+4),1) 500*us_bc(m_nl+4+1:end,type)] , 'LineStyle','non') ;
      set(h(1),'FaceColor',[1 1 1])  
    if type == 2  
      set(h(2),'FaceColor',[0.5,1,1])  % blue
    elseif type == 3  
       set(h(2),'FaceColor',[0.5,1,0.5]) % Green
    elseif type == 5  
       set(h(2),'FaceColor',[1, 0.75, 0.5])  % Yellow       
    end 
      
      L3 = h(2); 
%     hold on, 
%       h=area(ti,[ ( multp(nl+3+1:ns-2)-multps(nl+3+1:ns-2))' ( 2*multps(nl+3+1:ns-2))' ] );
%       set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
%  %      set(h(2),'FaceColor',[0.5 1 1])  % sky blue      
%       set(h(2),'FaceColor',[1 0.5 0.5])  % pink
%       set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value   
  hold on
        L1= plot(ti, multp(i,nl+4+1:ns),'LineStyle','-','Color','r', 'LineWidth',2.5);   
        L2= plot(ti, mean(multp(i,nl+4+1:ns)).*ones(size(ti,2),1),'LineStyle','--','Color','b', 'LineWidth',2.5);
    
        plot(ti, zeros(size(ti,2),1),'LineStyle','-','Color','k', 'LineWidth',1.5);
%     hold on
%     plot(ti, multp_ma(nl+3+1:ns-2)-multps_ma(nl+3+1:ns-2),'LineStyle','--','Color','r', 'LineWidth',1);
%     plot(ti, multp_ma(nl+3+1:ns-2)+multps_ma(nl+3+1:ns-2),'LineStyle','--','Color','r', 'LineWidth',1);
    hold off    
    xlim([1950 2020]);
    if j == 1 % y
        ylim([ 0  4]);
    elseif j == 2 % c
        ylim([ -1  2]); 
    elseif j == 3 % tax
        ylim([ -1  2]); 
    elseif j ==5 % inv
        ylim([ -1  2]);
    else
        ylim([ -0.5  1]);
    end    
   
%        ylim([ -5  5]);

    hold off
    title(['$', char(asshock(shock_seq(shock_gov))), ...
           '\uparrow\ \rightarrow\ ', ...  
           char(m_asvar(j)), '$ ; Sum of 0 - ', char(cha_imp(i))  ], 'interpreter', 'latex','FontSize',16)
  end
  
  if    type == 2
      legend([L1,L2,L3],{'Time-varying','Average','Recessions'})
  elseif type == 3
      legend('  ', 'Active Monetary Policy', 'Time-varying','Average')   
  elseif type == 4
      legend('  ', 'Active Fiscal Policy', 'Time-varying','Average')
  elseif type == 5
      legend('  ', 'Passive Fiscal Policy', 'Time-varying','Average')   
  end    
       
    est_date = datestr(date);   
    name = ['./Fig/Mutp_',num2str(nimp-1),'_',char(var_seq(j)),'_',est_date];
    saveas(h_(5000+10*j),name,'fig');      

    
end

% %  %  期間によるimpact fiscal multiplierの計算
% %  
%  vimp=[3 5 9 13 17 21];  % インパルスを合計する期間
% 
% cha_imp = {'0.5 year ahead';'1 year ahead';'2 years ahead';...
%            '3 years ahead'; '4 years ahead';'5 years ahead' };
% cha_imp_1 = {'2';'4';'8';'12'; '20' };
% 
% save_multp=zeros(size(vimp,2),nk,ns); 
% save_multps=zeros(size(vimp,2),nk,ns);
% multp=zeros(ns,1);
% multps=zeros(ns,1);
% 
% for i = 1: size(vimp,2) 
% % h_(5000+10*i) = figure(5000+10*i);
% 
%    for j = 1 : nk 
%      idg = (shock_seq(1)-1)*nk + 1;  
%     id = (shock_seq(1)-1)*nk + j;
%      mimp_g = reshape(mimpm(:, idg), nimp, ns)';     
%     mimp = reshape(mimpm(:, id), nimp, ns)';
%     
%      mimps_g = reshape(mimpms(:, idg), nimp, ns)';     
%     mimps = reshape(mimpms(:, id), nimp, ns)';
%     
%    
%     for k = nl+1:ns
%        
%        l= vimp(i); %nimp;  % インパルスを合計する期間
%         multp(k)=mimp(k,l)/mimp_g(k,1);         
%         multps(k)=mimps(k,l)/mimp_g(k,1);  
%       
%     end
%     
%      save_multp(i,j,:) = multp;
%      save_multps(i,j,:) = multps;
% 
%    end
%    
% end
% 
% for j=2:nk
%      h_(5000+10*j) = figure(5000+10*j);
%       set( h_(5000+10*j),'Position',[20,20,900,600]);
%       
%   for i = 1:size(vimp,2) 
%       
%       multp= squeeze(save_multp(:,j,:));
%       
%      subplot(nk/2,2, i);
% %       subplot(nk/3,3, j);
%      ti = 1953+(nl+4)/4:0.25:1953+(ns-2)/4;
% 
%         h=area(ti, [ -10*ones(239,1) 20*us_bc(6:(247-3),type)] , 'LineStyle','non') ;
%       set(h(1),'FaceColor',[1 1 1])  
%       set(h(2),'FaceColor',[0.5,1,1])
%         L3 = h(2);
% 
%       hold on
%         L1 = plot(ti, multp(i,nl+3+1:ns-2),'LineStyle','-','Color','r', 'LineWidth',2.5);
%         L2 = plot(ti, mean(multp(i,nl+3+1:ns-2)).*ones(size(ti,2),1),'LineStyle','--','Color','b', 'LineWidth',2.5);
%     
%         plot(ti, zeros(size(ti,2),1),'LineStyle','-','Color','k', 'LineWidth',1.5);
% 
%     hold off    
%     xlim([1950 2015]);
%     if j == 2
%         ylim([ -2  2]);
%     elseif j==3
%         ylim([ -1  1]);
%     elseif j == 4
%         ylim([ -0.5  3]);   
%     else
%         ylim([ -1  1]);   
%     end    
%        
%     hold off
%     title(['$', char(m_asvar(1)), ...
%            '\uparrow\ \rightarrow\ ', ...  
%            char(m_asvar(j)), '$ ; at ', char(cha_imp(i))  ], 'interpreter', 'latex','FontSize',12)
%   end
%   
%   if    type == 2
%       legend([L1,L2,L3],{'Time-varying','Average','Recessions'})
%   elseif type == 3
%       legend('  ', 'Active Monetary Policy', 'Time-varying','Average')   
%   elseif type == 4
%       legend('  ', 'Active Fiscal Policy', 'Time-varying','Average')
%   elseif type == 5
%       legend('  ', 'Passive Fiscal Policy', 'Time-varying','Average')   
%   end    
%        
%     est_date = datestr(date);   
%     name = ['./Fig/imp_Mutp_',num2str(nimp-1),'_',char(var_seq(j)),'_',est_date];
%     saveas(h_(5000+10*j),name,'fig');      
% 
%  end