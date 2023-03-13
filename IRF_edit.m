% IRF Results
%%--------------------------------------------------------%%
%%                     IRF_edit.m                       %%
%%--------------------------------------------------------%%



% clear all;

mimpm = load('tvpvar_imp.xls');
ns = 248; % # of time periods
nimp = size(mimpm, 1) / ns;

for i = 2;
    mimp = reshape(mimpm(:, i), nimp, ns)';
    save('tvpvar_imp_gy.xls', 'mimp',  '-ascii', '-tabs');

end;
for i = 3;
    mimp = reshape(mimpm(:, i), nimp, ns)';
    save('tvpvar_imp_gc.xls', 'mimp',  '-ascii', '-tabs');

end;
for i = 4;
    mimp = reshape(mimpm(:, i), nimp, ns)';
    save('tvpvar_imp_gpbd.xls', 'mimp',  '-ascii', '-tabs');

end;
for i = 5;
    mimp = reshape(mimpm(:, i), nimp, ns)';
    save('tvpvar_imp_gnex.xls', 'mimp',  '-ascii', '-tabs');

end;
for i = 6;
    mimp = reshape(mimpm(:, i), nimp, ns)';
    save('tvpvar_imp_grex.xls', 'mimp',  '-ascii', '-tabs');

end;