
%    jj =1;

  mimpm = [];   
    mimpms = []; 
    for i = 1:size(mimpm_temp,3)
        mimpm = [mimpm; mimpm_temp(:,:,i)/n_IRF(i) ] ;
        mimpms = [mimpms;  sqrt(mimpms_temp(:,:,i)/n_IRF(i)-(mimpm_temp(:,:,i)/n_IRF(i)).^2) ] ;
    end 
    
     %  財政ショックのサイズの基準化                     
    change_shock_size; 
    
    
    
     save(['./output/tvpvar_imps_' num2str(jj)  '.xls'], 'mimpms',  '-ascii', '-tabs');

     save( ['./output/tvpvar_imp_' num2str(jj) '.xls'], 'mimpm',  '-ascii', '-tabs');

    save_name_para = ['./output/para_save_' num2str(jj) '.mat'];
    save(save_name_para, 'nl','nimp', 'nsim', 'msampb', 'msampa', 'msamph','n_IRF' );  