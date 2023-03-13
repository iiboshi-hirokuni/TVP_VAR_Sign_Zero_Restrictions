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
% filename=('./save_output/'); 

mimpm = load([ filename '/tvpvar_imp_',char(policy_type),'-' ,char(svar_type),'.xls']);     % posterior means of IRF
mimpms = load( [ filename '/tvpvar_imps_',char(policy_type),'-' ,char(svar_type),'.xls' ] );   % Standard Deviation of IRF

m_asvar = {'y';'tax';'g';'int';'\pi'};                       
%%  インパルス応答の作成

% policy = 'Tax'
 policy = 'FP'

 j = 1;   % 出力したい変数; 1:'Y';2'C';3'Tax'; 4'Gov';5 'Inv';6 'R';7 'Pi'

global m_ns m_nk m_nl m_asvar m_vpd policy_type;

ns = m_ns;  
nk = m_nk; 
nl = m_nl; 
asvar = {'y';'tax';'g';'int';'\pi'};   

var_seq = { 'Y';'Tax'; 'Gov'; 'R'; 'pi' };

switch policy
    case 'FP'
        switch svar_type
          case 'Mountford_Uhlig'   
            var_Gov = 4;
            shock_gov = 3;
            shock_seq =[1 2 3 4];  % 1:BC shock, 2:MP shock, 3:FP shock, 4:TAX 
            asshock = {'BC';'MP';'FP';'TAX';' ';'real';'Demand';'non'};
            graph_start= 1;
            graph_end  = 4;     
               idx = 2;    
             shock_s =3;
             shock_e =3;  
          case 'Caldala_Kamp' 
            var_Gov = 3;
            shock_gov = 1;
            shock_seq =[1 2 3 4];  % 1:FP shock 
            asshock = {'gov';'Tax';' ';'real';'Demand';'non'};
             graph_start= 1;
             graph_end  = 4;     
               idx = 2;    
             shock_s =1;
             shock_e =1; 
         case 'Caldala_Kamps_tax' 
            var_Gov = 2;
            shock_gov = 3;
            shock_seq =[1 2 3 4];  % 1:FP shock 
            asshock = {'gov';' ';'Tax Cut';'Demand';'non'};
             graph_start= 1;
             graph_end  = 4;     
               idx = 2;    
             shock_s =3;
             shock_e =3;   
             
             mimpm = -1*mimpm;
         end    
      
    case 'Tax'
          var_Gov = 2;
          shock_gov = 4;
end          
 

nimp = size(mimpm, 1) / m_ns;
mline = [1 0 0;  0 0 1; 0 .5 0;  0 .7 .7];
vline = {':.', '-', '--', ':'};
% nline = size(vt, 2);

switch  policy_type
   case 'Choleski'
    asshock = {'BC';'MP';'FP';'TAX';' ';'real';'Demand';'non'};
    graph_start= 1;
    graph_end  = 4;
    idx = 2;     
     shock_seq =[1 2 3 4];  % 1:BC shock, 2:MP shock, 3:FP shock, 4:TAX  
end

time_title={'1960Q1' '1970Q1' '1980Q1' '1990Q1' '2000Q1' '2010Q1' }'; 
time_seq  =[37 77 117 157 197 237]; 

ti = 1952:0.25:1952+(ns-nl-1)/4;
horizon =[1 2 4 8 12 20];
hz = 0:1:(nimp-1); 


 %  インパルス応答の作成
 
for j = 1:5 %出力したい変数; 'Y';'Tax'; 'Gov';'R'; 'Pi' 
 for i = shock_s: shock_e  % shock, 1:Gov, 3:Tax  
     
        
    h_(3000+i+10*j) =figure(3000+i+10*j);
       set( h_(3000+i+10*j),'Position',[20,20,900,600]);
    
      id = (shock_seq(i)-1)*nk + j;
      mimp = reshape(mimpm(:, id), nimp, ns)';
      mimps = reshape(mimpms(:, id), nimp, ns)';   
     
     y_max= 1.2*max(max(mimp+mimps)');
     y_min= 1.2*min(min(mimp-mimps)');
   
    
    for k = 1:6
     subplot(4,3,k)   
%      id = (shock_seq(i)-1)*nk + j;
%      mimp = reshape(mimpm(:, id), nimp, ns)';
%      mimps = reshape(mimpms(:, id), nimp, ns)'   
    
      h=area(hz,[ (mimp(time_seq(k)+1,:)-mimps(time_seq(k)+1,:))' (2*mimps(time_seq(k)+1,:))' ] );
      set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
      set(h(2),'FaceColor',[0.5 1 1])      
      set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
      hold on
      plot(hz,mimp(time_seq(k)+1,:),'LineStyle','-','Color','b',...
        'LineWidth',2.5),  
       plot(hz,0*mimp(time_seq(k)+1,:),'LineStyle','-','Color','k',...
        'LineWidth',0.5),  
      hold off
%       title([ time_title(k), char(m_asvar(j)),  char(asshock(i)) ])
       title([ char(var_seq(j)), ' to ' , char(asshock(i)), ' shock, ' , char(time_title(k)) ],...
           'FontSize',11,'FontWeight','bold');
%       xlim([0 nimp-1 ])
      xlim([0 20 ])      
      ylim([y_min y_max] )
    end 
    
    for l = 1:6
     subplot(4,3,6+l)
      
      hold on
      h=area(ti,[ (mimp(nl+1:end,horizon(l)+1)-mimps(nl+1:end,horizon(l)+1)) ( 2*mimps(nl+1:end,horizon(l)+1)) ] );
      set(h(1),'FaceColor',[1 1 1])        % [.5 0 0])
%       set(h(2),'FaceColor',[0.5 1 1])  % sky blue      
      set(h(2),'FaceColor',[1 0.5 0.5])  % pink
      set(h,'LineStyle','none','LineWidth',0.5) % Set all to same value
     
     plot (ti,mimp(nl+1:end,horizon(l)+1), 'LineStyle','-','Color','r',...
        'LineWidth',2.5),  
     plot (ti,0*mimp(nl+1:end,horizon(l)+1), 'LineStyle','-','Color','k',...
        'LineWidth',0.5),  

      hold off
     xlim( [1950 2020  ]);
      ylim([y_min y_max] )

     title( [ num2str(horizon(l)), ' period ahead'],'FontSize',11,'FontWeight','bold');
    end   
    
    est_date = datestr(date);   
    name = ['./Fig/imp_',num2str(nimp-1),'_',char(asshock(i)),'_',char(var_seq(j)),'_',est_date];
    saveas(h_(3000+i+10*j),name,'fig'); 
     
 end
  
 end
 
%  
%  %  期間による累積の財政乗数の計算
%  
%  vimp=[3 5 9 13 21];  % インパルスを合計する期間
% 
% cha_imp = {'2 \ periods';'4 \ periods';'8 \ periods';'12 \ periods'; '20 \ periods' };
% 
% for i = 1: size(vimp,2) 
%   h_(3000+i+20*j) = figure(3000+i+20*j);
%       set( h_(3000+i+20*j),'Position',[20,20,900,600]);
% 
% 
%    for j = 1 : nk 
%      idg = (shock_seq(shock_gov)-1)*nk + var_Gov;  
%     id = (shock_seq(shock_gov)-1)*nk + j;
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
%         multp(k)=sum(mimp(k,1:l),2)/sum(mimp_g(k,1:l),2);        
%          multps(k)=sum(mimps(k,1:l),2)/sum(mimp_g(k,1:l),2);       
% %        multps(k)=mean(mimps(k,1:l),2) ;  
%       
%     end
%     
%        
%      subplot(ceil(nk/2),2, j);
% %       subplot(nk/3,3, j);
%      ti = 1953+(nl+4)/4:0.25:1953+(ns-2)/4;
%      
%     hold on
%            plot(ti, multp(nl+3+1:ns-2),'LineStyle','-','Color','r', 'LineWidth',2.5);
%         
%     hold on
%         plot(ti, mean(multp(nl+3+1:ns-2)).*ones(size(ti,2),1),'LineStyle','--','Color','b', 'LineWidth',2.5);    
%         plot(ti, zeros(size(ti,2),1),'LineStyle','-','Color','k', 'LineWidth',1.5);
% 
%     hold off    
%     xlim([1960 2020]);
%     if j == 2
%         ylim([ -3  4]);
%     end    
%        
%     hold off
%     title(['$', char(var_seq(var_Gov)), ...
%            '\uparrow\ \rightarrow\ ', ...
%            char(var_seq(j)), ' ; 0 \ to \ ', char(cha_imp(i)),  '$'], 'interpreter', 'latex','FontSize',12)
%            
% 
%   end
%   
% end
 