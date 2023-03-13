%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  [] = drawimp(vt, fldraw)
%%
%%  "drawimp" draws time-varying impulse response
%%
%%  [input]
%%   (fldraw = 1)
%%     vt:   m*1 vector of horizons to draw impulse
%%   (fldraw = 0)
%%     vt:   m*1 vector of time points to draw impulse
%%


% discount factor
beta = 0.995; %0.995; 


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
           % 6: Political Regime 
           
% load('para_save.mat');

us_bc=csvread('./data/US_BC.csv',1);



nimp = size(mimpm, 1) / m_ns;
mline = [1 0 0;  0 0 1; 0 .5 0;  0 .7 .7];
vline = {':.', '-', '--', ':'};
% nline = size(vt, 2);

switch  policy_type
  case 'unanticipated_fiscal_policy' 
     asshock = {'FP';'BC';'MP';'real';'Demand';'non'};
%      graph_start= 1;
%      graph_end  = 3;
      idx = 1; 
     shock_seq =[1 2 3];  % 1:FP shock, 2:BC shock, 6:MP shock 
  case 'anticipated_fiscal_policy'  
    asshock = {'';'';'FP';'BC';'MP';'real';'Demand';'non'};
%     graph_start= 3;
%     graph_end  = 5;
     idx = 3;
 case 'bench_mark'
    asshock = {'FP';'BC';'MP';' ';'real';'Demand';'non'};
     graph_start= 1;
     graph_end  = 3;
     idx = 2;    
    shock_seq =[3 4 2];  % 1:FP shock, 2:BC shock, 6:MP shock 
end

time_title={'1960Q1' '1970Q1' '1980Q1' '1990Q1' '2000Q1' '2010Q1' }'; 
time_seq  =[33 73 113 153 193 233]; 

ti = 1953:0.25:1953+(ns-nl-1)/4;
horizon =[1 4 8 12 16 20];
hz = 0:1:(nimp-1); 

var_seq = { 'Gov'; 'Y'; 'Cons'; 'Debt'; 'Price'; 'R' };

 
 
 %  期間による累積の財政乗数の計算
 
 vimp=[3 5 13  21];  % インパルスを合計する期間
 nfig = length(vimp) 

cha_imp = {'0.5 yrs';'1 yrs';...
           '3 yrs'; '5 yrs' };
cha_imp_1 = {'2';'4';'8';'12'; '20' };

save_multp=zeros(size(vimp,2),nk,ns); 
save_multps=zeros(size(vimp,2),nk,ns);
multp=zeros(ns,1);
multps=zeros(ns,1);


for kk = 1:3
    
  if kk == 1  
      filename=filename1;  % ('../Waseda_results/output_with_zero_V1/'); 
  elseif kk == 2
       filename=filename2; % ('../Waseda_results/output_only_sign_V1/'); 
  else
       filename=filename3;
  end
       mimpm = load([ filename '/tvpvar_imp.xls']);     % posterior means of IRF

for i = 1: size(vimp,2) 

   for j = 1 : nk 
     idg = (shock_seq(1)-1)*nk + 1;  
    id = (shock_seq(1)-1)*nk + j;
     mimp_g = reshape(mimpm(:, idg), nimp, ns)';     
    mimp = reshape(mimpm(:, id), nimp, ns)';
    
     mimps_g = reshape(mimpms(:, idg), nimp, ns)';     
    mimps = reshape(mimpms(:, id), nimp, ns)';
    
   
     for k = nl+1:ns
       
       l= vimp(i); %nimp;  % インパルスを合計する期間
%         multp(k)=sum(mimp(k,1:l),2)/sum(mimp_g(k,1:l),2);   
        %% 乗数の現在割引価値 
         multp(k)=sum_pv(mimp(k,1:l),beta)/sum_pv(mimp_g(k,1:l),beta);  
           
         multps(k)=sum(mimps(k,1:l),2)/sum(mimp_g(k,1:l),2);         
      
     end
    
     save_multp(i,j,:) = multp;
     save_multps(i,j,:) = multps;
    
   end      
 
end
  
   if kk == 1 
      save_multp1= save_multp;
   elseif kk== 2       
      save_multp2= save_multp;
   else
      save_multp3= save_multp;
   end 
 
end


for j=2:nk
     h_(6000+10*j) = figure(6000+10*j);
      set( h_(6000+10*j),'Position',[20,20,900,600]);
      
  for i = 1:size(vimp,2) 
      
      multp1= squeeze(save_multp1(:,j,:));
      multp2= squeeze(save_multp2(:,j,:));
      multp3= squeeze(save_multp3(:,j,:));
      
     subplot(nfig/2,2, i);
%       subplot(nk/3,3, j);
     ti = 1953+(nl+4)/4:0.25:1953+(ns-2)/4;
     
     h=area(ti, [ -20*ones(239,1) 50*us_bc(6:(247-3),type)] , 'LineStyle','non') ;
      set(h(1),'FaceColor',[1 1 1])  
    if type == 2  
      set(h(2),'FaceColor',[0.5,1,1])  % blue
    elseif type == 3  
       set(h(2),'FaceColor',[0.5,1,0.5]) % Green
    elseif type == 5  
        set(h(2),'FaceColor',[1, 0.75, 0.5])  % Yellow            
    elseif type == 6  
         set(h(2),'FaceColor',[1,0.5,1])  % blue
    end 
      
      L4 = h(2); 
%     hold on, 
%       h=area(ti,[ ( multp(nl+3+1:ns-2)-multps(nl+3+1:ns-2))' ( 2*multps(nl+3+1:ns-2))' ] );
%       set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
%  %      set(h(2),'FaceColor',[0.5 1 1])  % sky blue      
%       set(h(2),'FaceColor',[1 0.5 0.5])  % pink
%       set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value   
  hold on
        L1= plot(ti, multp1(i,nl+3+1:ns-2),'LineStyle','-','Color','r', 'LineWidth',2.5);        
        L2= plot(ti, multp2(i,nl+3+1:ns-2),'LineStyle','--','Color','b', 'LineWidth',2.5); 
        L3= plot(ti, multp3(i,nl+3+1:ns-2),'LineStyle',':','Color','k', 'LineWidth',2.5); 
%         L2= plot(ti, mean(multp(i,nl+3+1:ns-2)).*ones(size(ti,2),1),'LineStyle','--','Color','b', 'LineWidth',2.5);
    
        plot(ti, zeros(size(ti,2),1),'LineStyle','-','Color','k', 'LineWidth',1.5);
%     hold on
%     plot(ti, multp_ma(nl+3+1:ns-2)-multps_ma(nl+3+1:ns-2),'LineStyle','--','Color','r', 'LineWidth',1);
%     plot(ti, multp_ma(nl+3+1:ns-2)+multps_ma(nl+3+1:ns-2),'LineStyle','--','Color','r', 'LineWidth',1);
    hold off    
    xlim([1950 2015]);
    if j == 2 % y
        ylim([ -1  2]);
    elseif j == 3 % c
        ylim([ -2  2]); 
    elseif j == 4 % debt
        ylim([ 0  5]); 
    elseif j ==5 % inf
        ylim([ -2  2]);
    else
        ylim([ -1  1]);
    end    
       
    hold off
    title(['$', char(m_asvar(1)), ...
           '\uparrow\ \rightarrow\ ', ...  
           char(m_asvar(j)), '$ ; Cusum b/w 0 \& ', char(cha_imp(i))  ], 'interpreter', 'latex','FontSize',12)
  end
  
  if    type == 2
%       legend([L1,L2,L3],{'w/ Zero','only Sign','Recessions'})
       legend([L1,L2,L3,L4],{'ALL shocks','FP and BC shocks', 'only FP shocks','Recessions'},'FontSize',11)
  elseif type == 3
      legend('  ', 'Active Monetary Policy', 'Time-varying','Average')   
  elseif type == 4
      legend('  ', 'Active Fiscal Policy', 'Time-varying','Average')
  elseif type == 5
      legend('  ', 'Passive Fiscal Policy', 'Time-varying','Average')   
  elseif type == 6
%       legend([L1,L2,L3],{'w/ Zero','only Sign','Recessions'})
       legend([L1,L2,L3,L4],{'ALL shocks','FP and BC shocks', 'only FP shocks','Republic'})  
  end    
       
    est_date = datestr(date);   
    name = ['./Fig/Comp_Mutp_',num2str(nimp-1),'_',char(var_seq(j)),'_',est_date];
    saveas(h_(6000+10*j),name,'fig');      

    
end

