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

% filename1=('./output/');
% filename2=('./output/');


set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');


% discount factor
avri = 0.04314;        %average of nominal interest rate
beta = 1/(1+avri); 

m_asvar = {'real GDP';'Tax Revenue ';'Gov. Expenditure';'Interest Rate';'Inflation'};
 
   %   Y/G            Tax/G  G/G   int/G  ing/g
% ratio_g =[ 4.77476,    1/1,   1/1,   1,     1    ];


%%  �C���p���X�����̍쐬

 j = 1;   % �o�͂������ϐ�; 1:Gov, 2:Y, 3:cons, 4:debt, 5:pi, 6:int 

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



mline = [1 0 0;  0 0 1; 0 .5 0;  0 .7 .7];
vline = {':.', '-', '--', ':'};
% nline = size(vt, 2);

switch  policy_type
%   case 'unanticipated_fiscal_policy' 
%      asshock = {'FP';'BC';'MP';'real';'Demand';'non'};
% %      graph_start= 1;
% %      graph_end  = 3;
%       idx = 1; 
%      shock_seq =[1 2 3];  % 1:FP shock, 2:BC shock, 6:MP shock 
%   case 'anticipated_fiscal_policy'  
%     asshock = {'';'';'FP';'BC';'MP';'real';'Demand';'non'};
% %     graph_start= 3;
% %     graph_end  = 5;
%      idx = 3;
 case 'bench_mark'
    asshock = {'Gov';'BC';'MP';' ';'real';'Demand';'non'};
     graph_start= 1;
     graph_end  = 3;
     idx = 2;    
    shock_seq =[3 4 2];  % 1:FP shock, 2:BC shock, 6:MP shock 
end

time_title={'1960Q1' '1970Q1' '1980Q1' '1990Q1' '2000Q1' '2010Q1' }'; 
time_seq  =[33 73 113 153 193 233]; 

ti = 1953:0.25:1953+(ns-nl-1)/4;
horizon =[1 4 8 12 16 20];


var_seq = { 'Y';'Tax'; 'Gov'; 'R'; 'Pi' };

 
 
 %  ���Ԃɂ��ݐς̍����搔�̌v�Z
 
 vimp=[ 5   21  33];  % �C���p���X�����v�������
 nfig= length(vimp);

cha_imp = { '1 yrs','5 yrs','8 yrs' };
cha_imp_1 = {'2';'4';'8';'12'; '20' };

save_multp=zeros(size(vimp,2),nk,ns); 
save_multps=zeros(size(vimp,2),nk,ns);
multp=zeros(ns,1);
multps=zeros(ns,1);


for kk = 1:2
    
  if kk == 1  
      filename=filename1;  % ('../Waseda_results/output_with_zero_V1/'); 
       var_Gov = 3;
        shock_seq =[1 2 3];
        shock_gov =1; 
        svar_type = 'Caldala_Kamp'
               %   Y/G         Tax/G  G/G  int/G  inf/g
            ratio_g =[ 4.83099,  1/1,   1/1,  1,     1    ];
  elseif kk == 2
       filename=filename2; % ('../Waseda_results/output_only_sign_V1/'); 
         var_Gov = 2; % tax
         shock_gov = 3;           
         shock_seq =[1 2 3 4];  % 1:FP shock 
         svar_type = 'Caldala_Kamps_tax'
         
          ratio_g =[ 5.670985,  1/1,   1/1,  1,     1    ];
  elseif kk == 3
       filename=filename3; % ('../Waseda_results/output_only_sign_V1/');    
  else
       filename=filename4;
  end

       mimpm = load([ filename '/tvpvar_imp_',char(policy_type),'-' ,char(svar_type),'.xls']);     % posterior means of IRF
       mimpms = load([ filename '/tvpvar_imps_',char(policy_type),'-' ,char(svar_type),'.xls']);     % posterior means of IRF
       nimp = size(mimpm, 1) / m_ns;
       hz = 0:1:(nimp-1); 
       
 for i = 1: size(vimp,2) 

   for j = 1 : nk 
    idg = (shock_seq(shock_gov)-1)*nk + var_Gov;  
    id = (shock_seq(shock_gov)-1)*nk + j;
     mimp_g = reshape(mimpm(:, idg), nimp, ns)';     
    mimp = reshape(mimpm(:, id), nimp, ns)';
    
     mimps_g = reshape(mimpms(:, idg), nimp, ns)';     
    mimps = reshape(mimpms(:, id), nimp, ns)';
    
   
     for k = nl+1:ns
       
       l= vimp(i); %nimp;  % �C���p���X�����v�������
%         multp(k)=sum(mimp(k,1:l),2)/sum(mimp_g(k,1:l),2);   
        %% �搔�̌��݊������l 
       if kk == 2 
         multp(k)=(-1)*sum_pv(mimp(k,1:l),beta)/sum_pv(mimp_g(k,1:l),beta)*ratio_g(j); 
       else
         multp(k)=sum_pv(mimp(k,1:l),beta)/sum_pv(mimp_g(k,1:l),beta)*ratio_g(j); 
       end    
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
   elseif kk== 3       
      save_multp3= save_multp;
   elseif kk== 4
      save_multp4= save_multp;
   end 
 
end


for j=1:3 %nk
     h_(6000+10*j) = figure(6000+10*j);
      set( h_(6000+10*j),'Position',[20,20,900,600]);
      
  for i = 1:size(vimp,2) 
      
      multp1= squeeze(save_multp1(:,j,:));
      multp2= squeeze(save_multp2(:,j,:));
      
     subplot(size(vimp,2),1, i);
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
    elseif type == 6  
         set(h(2),'FaceColor',[1,0.5,1])  % blue
    end 
      
      L4 = h(2); 
 
  hold on
        L1= plot(ti, multp1(i,nl+4+1:ns),'LineStyle','-','Color','r', 'LineWidth',1.5); 
        L3= plot(ti, mean(multp1(i,nl+4+1:ns),2).*ones(size(multp1(i,nl+4+1:ns),2),1),'LineStyle',':','Color','r', 'LineWidth',3); 
        L2= plot(ti, multp2(i,nl+4+1:ns),'LineStyle','-','Color','b', 'LineWidth',1.5); 
        L4= plot(ti, mean(multp2(i,nl+4+1:ns),2).*ones(size(multp1(i,nl+4+1:ns),2),1),'LineStyle','--','Color','b', 'LineWidth',3); 
         plot(ti, zeros(size(ti,2),1),'LineStyle','-','Color','k', 'LineWidth',1.5);

    hold off    
    xlim([1950 2020]);
    if j == 1 % y
        ylim([ -1  3]);
    elseif j == 2 % gov
        ylim([ -1.5  2]); 
    elseif j == 3 % tax
        ylim([ -3  2]); 
    elseif j ==5 % pi
        ylim([ -4  4]);
    else
        ylim([ -4  4]);
    end  
  
%       ylim([ -2  3 ]);
       
    hold off
%     title(['$', char(asshock(1)), ...
     title([' Fiscal Policy',  '$ \rightarrow\ $', ...  
           ' Cumulative Effect of ', char(m_asvar(j)), ',  For ', char(cha_imp(i))  ], 'interpreter', 'latex','FontSize',14)
%            char(m_asvar(j)), '$ ; Cusum b/w 0 \& ', char(cha_imp(i))  ], 'interpreter', 'latex','FontSize',14)
             
  end
  
  if    type == 2
%          legend([L1,L2,L3,L4],{'Benchmark','Model 2','Model 3','Model 4'})
          legend([L1,L3,L2,L4],{'Gov. Spend. Rule','Mean of GSR','Tax Cut Rule','Mean of TCR'},'FontSize',11)
%       legend([L1,L2,L3],{'w/ Zero','only Sign','Recessions'})
%        legend([L1,L2,L3,L4],{'ALL shocks','FP and BC shocks', 'only FP shocks','Recessions'})
  
  end    
       
    est_date = datestr(date);   
    name = ['./Fig/Comp_Mutp_',char(var_seq(j)),'_',est_date];
    saveas(h_(6000+10*j),name,'fig');      

    
end

