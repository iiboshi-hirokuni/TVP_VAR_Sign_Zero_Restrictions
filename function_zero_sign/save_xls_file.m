

est_date = datestr(date)

mimpm =load('./output/tvpvar_imp.xls');


file_name = ['./output/tvpvar_imp_',num2str(jj) ,'_' ,num2str(nsim), '_', est_date,'_', policy_type, svar_type, '.xls']

save(file_name, 'mimpm', '-ascii', '-tabs');

file_name_s = ['./output/tvpvar_imps_',num2str(jj) ,'_' ,num2str(nsim), '_', est_date,'_', policy_type,  '.xls']
save(file_name_s, 'mimpms', '-ascii', '-tabs');