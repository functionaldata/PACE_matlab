function[ev_trunc, ev_fit] = myEigenvalues(dat_all, t_all, no_eig, p1)
% dat_all:  a 1 by N list of function values, where N is the number of
% subjects.
% t_all:    a 1 by N list of the corresponding time points.
% no_eig:   number of eigenvalues to obtain for both methods.
% p1:       settings for FPCA

tmp = FPCA(dat_all, t_all, p1);
rcov = tmp{37};
phi = tmp{4};
ev_trunc = tmp{3}(1:no_eig);
out1 = tmp{20};

y = rcov.cxxn;
X = getDesignMat(rcov, phi, out1, no_eig);
matlab_ver = version('-release');
% Matlab version problem
if (strcmp(matlab_ver, '2013b') || str2num(matlab_ver(1:4)) > 2013)
    model = fitlm(X, y, 'Intercept', false);
else
    model = LinearModel.fit(X, y, 'Intercept', false);
end

ev_fit = double(model.Coefficients(:, 1))';
% plot(tmp{4}(:, 2))

end

function[X] = getDesignMat(rcov, phi, out1, no_eig)

lam1 = interp1(out1, phi(:, 1:no_eig), rcov.tpairn(1, :));
lam2 = interp1(out1, phi(:, 1:no_eig), rcov.tpairn(2, :));
X = lam1 .* lam2;

end