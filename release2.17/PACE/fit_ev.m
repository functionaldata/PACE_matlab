function [ ev_fit ] = fit_ev( rcov, phi, out21, no_eig )
%fit_ev Use least-squares fit method to get the eigenvalues.
%   Input the raw covariance matrix rcov, eigenfunctions phi, evaluation
%   time points out21, and number of eigenvalues to fit no_eig

y = rcov.cxxn;
X = getDesignMat(rcov, phi, out21, no_eig);
matlab_ver = version('-release');
% Matlab version problem
if (strcmp(matlab_ver, '2013b') || str2num(matlab_ver(1:4)) > 2013)
    model = fitlm(X, y, 'Intercept', false);
else
    model = LinearModel.fit(X, y, 'Intercept', false);
end

ev_fit = model.Coefficients.Estimate;

end


function[X] = getDesignMat(rcov, phi, out1, no_eig)

lam1 = interp1(out1, phi(:, 1:no_eig), rcov.tpairn(1, :));
lam2 = interp1(out1, phi(:, 1:no_eig), rcov.tpairn(2, :));
X = lam1 .* lam2;

end