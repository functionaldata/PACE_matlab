ms = table2array(readtable('ev_trunc_mild_std.txt', 'Delimiter', '\t', 'ReadVariableNames', false));
mt = table2array(readtable('ev_trunc_mild_true_mean.txt', 'Delimiter', '\t', 'ReadVariableNames', false));
vs = table2array(readtable('ev_trunc_valley_std.txt', 'Delimiter', '\t', 'ReadVariableNames', false));
vt = table2array(readtable('ev_trunc_valley_true_mean.txt', 'Delimiter', '\t', 'ReadVariableNames', false));

lambda_true = 1 ./ (1:5).^2;
no_eig = 5;
mean_mis = mean(vs);
std_mis = std(vs) / sqrt(size(ms, 1));
out_cell = cell(1, 5);
for (j = 1:no_eig)
    out_cell{j} = sprintf('%6.3f (%4.3f)',  mean_mis(j) - lambda_true(j), std_mis(j));
end
out_cell

M = size(ms, 1);
mean_mis = mean((1 ./ vt - ... 
    ones(size(ms, 1), 1) * (1 ./ lambda_true(1:no_eig))) .^ 2)
std_mis = std((1 ./ vt - ... 
    ones(size(ms, 1), 1) * (1 ./ lambda_true(1:no_eig))) .^ 2) / sqrt(M)    
for (j = 1:no_eig)
    out_cell{j} = sprintf('%6.3f',  mean_mis(j));
end
out_cell
