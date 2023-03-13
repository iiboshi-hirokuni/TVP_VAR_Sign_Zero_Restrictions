function [pv_cusum] = sum_pv(series, beta)

num =size(series,2);

total = series(1);

for i = 2:num
    total = total + beta^i*series(i);
end    

pv_cusum = total;