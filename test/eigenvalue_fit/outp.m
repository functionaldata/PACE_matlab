no_eig = size(ev_trunc0, 2);
M = size(ev_trunc0, 1);
bias = nan * ones(4, no_eig);
mse = nan * ones(4, no_eig);
mrd = mse;
std_mat = mse;

means = lambda_true(1:no_eig);

bias(1, :) = mean(ev_trunc1) - means;
bias(2, :) = mean(ev_trunc0) - means;
bias(3, :) = mean(ev_fit1) - means;
bias(4, :) = mean(ev_fit0) - means;

std_mat(1, :) = std(ev_trunc1) / sqrt(M);
std_mat(2, :) = std(ev_trunc0) / sqrt(M);
std_mat(3, :) = std(ev_fit1) / sqrt(M);
std_mat(4, :) = std(ev_fit0) / sqrt(M);

names = {'trunc1', 'trunc0', 'fit1', 'fit0'};

mse(1, :) = mean((ev_trunc1 - ones(M, 1) * means).^2);
mse(2, :) = mean((ev_trunc0 - ones(M, 1) * means).^2);
mse(3, :) = mean((ev_fit1 - ones(M, 1) * means).^2);
mse(4, :) = mean((ev_fit0 - ones(M, 1) * means).^2);

std_mse(1, :) = std((ev_trunc1 - ones(M, 1) * means).^2) / sqrt(M);
std_mse(2, :) = std((ev_trunc0 - ones(M, 1) * means).^2) / sqrt(M);
std_mse(3, :) = std((ev_fit1 - ones(M, 1) * means).^2) / sqrt(M);
std_mse(4, :) = std((ev_fit0 - ones(M, 1) * means).^2) / sqrt(M);

mrd(1, :) = mean(abs(relDiff(ev_trunc1, means)));
mrd(2, :) = mean(abs(relDiff(ev_trunc0, means)));
mrd(3, :) = mean(abs(relDiff(ev_fit1, means)));
mrd(4, :) = mean(abs(relDiff(ev_fit0, means)));

out_cell0 = cell(1, 6);
out_cell0{1, 1} = 'mean';
for (j = 2:size(out_cell0, 2))
    out_cell0{1, j} = sprintf('%6.4f', means(j - 1));
end

out_cell1 = cell(5, 6);
out_cell1{1, 1} = 'bias';
for (i = 2:size(out_cell1, 1))
    out_cell1{i, 1} = names{i - 1};
    for (j = 2:size(out_cell1, 2))
        out_cell1{i, j} = sprintf('%6.4f (%4.4f)', bias(i - 1, j - 1), std_mat(i - 1, j - 1));
    end
end

out_cell2 = cell(5, 6);
out_cell2{1, 1} = 'MSE';
for (i = 2:size(out_cell2, 1))
    out_cell2{i, 1} = names{i - 1};
    for (j = 2:size(out_cell2, 2))
        out_cell2{i, j} = sprintf('%6.4f', mse(i - 1, j - 1));
    end
end

out_cell3 = cell(5, 6);
out_cell3{1, 1} = 'MRD';
for (i = 2:size(out_cell3, 1))
    out_cell3{i, 1} = names{i - 1};
    for (j = 2:size(out_cell3, 2))
        out_cell3{i, j} = sprintf('%6.4f', mrd(i - 1, j - 1));
    end
end

out_cell = [out_cell0; out_cell1; out_cell2; out_cell3];
xlswrite('tmp.xls', out_cell)

