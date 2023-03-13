%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');




% discount factor
beta = 0.995; 

        %   Y/G      Tax/G  G/G  int/G  inf/g
ratio_t =[ 5.670985,  1/1,   1/1,  1,     1    ];

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

switch  policy_type  
 case 'bench_mark'
     switch svar_type        
         case 'Caldala_Kamps_tax' 
           var_seq = { 'Y';'Tax'; 'Gov'; 'R'; 'Pi' };
            var_tax = 2; % tax
            var_gov = 2; % tax
            shock_tax = 3;           
            shock_seq =[1 2 3 4];  % 1:FP shock 
            asshock = {'gov';' ';'Tax';'Demand';'non'};
             graph_start= 1;
             graph_end  = 4;     
               idx = 2;    
             shock_s =3;
             shock_e =3;   
         end    
end

time_title={'1960Q1' '1970Q1' '1980Q1' '1990Q1' '2000Q1' '2010Q1' }'; 
time_seq  =[33 73 113 153 193 233]; 

ti = 1952:0.25:1952+(ns-nl-1)/4;
horizon =[1 4 8 12 16 20];
hz = 0:1:(nimp-1); 

  
 
 %  期間による累積の財政乗数の計算
 
  vimp=[3  5  13 21];  % インパルスを合計する期間

cha_imp = {'0.5 years';'1 years';...
           '3 years'; '5 years' };
cha_imp_1 = {'2';'4';'8';'12'; '20' };;

save_multp=zeros(size(vimp,2),nk,ns); 
save_multps=zeros(size(vimp,2),nk,ns);
multp=zeros(ns,1);
multps=zeros(ns,1);

for i = 1: size(vimp,2) 
% h_(4000+10*i) = figure(4000+10*i);

   for j = 1 : nk 
     idg = (shock_seq(shock_tax)-1)*nk +  var_tax;  
    id = (shock_seq(shock_tax)-1)*nk + j;
     mimp_g = reshape(mimpm(:, idg), nimp, ns)';     
    mimp = reshape(mimpm(:, id), nimp, ns)';
    
     mimps_g = reshape(mimpms(:, idg), nimp, ns)';      
    mimps = reshape(mimpms(:, id), nimp, ns)';
    
   
    for k = nl+1:ns
       
       l= vimp(i); %nimp;  % インパルスを合計する期間
%         multp(k)=sum(mimp(k,1:l),2)/sum(mimp_g(k,1:l),2);       
       
        %% 乗数の現在割引価値 
         multp(k)= ((-1) *sum_pv(mimp(k,1:l),beta)/sum_pv(mimp_g(k,1:l),beta))*ratio_t(j);  

         multps(k)=sum(mimps(k,1:l),2)/sum(mimp_g(k,1:l),2);  
        
    end
    
     save_multp(i,j,:) = multp;
     save_multps(i,j,:) = multps;
    

   end
   
end

for j=1:nk
     h_(4000+10*j) = figure(4000+10*j);
      set( h_(4000+10*j),'Position',[20,20,900,600]);
      
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

  hold on
        L1= plot(ti, multp(i,nl+4+1:ns),'LineStyle','-','Color','r', 'LineWidth',2.5);   
        L2= plot(ti, mean(multp(i,nl+4+1:ns)).*ones(size(ti,2),1),'LineStyle','--','Color','b', 'LineWidth',2.5);
    
        plot(ti, zeros(size(ti,2),1),'LineStyle','-','Color','k', 'LineWidth',1.5);

    hold off    
    xlim([1950 2020]);
    if j == 1 % y
        ylim([ -4  4]);
    elseif j == 2 % c
        ylim([ -2  2]); 
    elseif j == 4 % gov
        ylim([ -3  1]); 
    elseif j ==5 % inv
        ylim([ 0  8]);
    else
        ylim([ -0.5  1]);
    end    
       
    hold off
    title(['$', char(asshock(shock_tax)), ...
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
    name = ['./Fig/Mutp_',num2str(nimp-1),'_',char(var_seq(j)),'_',char(svar_type),est_date];    
    saveas(h_(4000+10*j),name,'fig');      

    
end
