function[dat_thin, t_thin] = thinObs(dat_all, t_all, unif_on)
% The number of observations N is chosen discrete uniformly on unif_on.
% Then for each individuals N random samples of observations will be
% retained.

if (length(unif_on) == 1)
    unif_on = [unif_on, unif_on];
end

n = length(dat_all);
dat_thin = cell(1, n);
t_thin = dat_thin;
for (i = 1:n)
    n_orig = length(dat_all{i});
    n_thin = min(randsample(unif_on, 1), n_orig);
    ind = sort(randsample(n_orig, n_thin));
    dat_thin{i} = dat_all{i}(ind);
    t_thin{i} = t_all{i}(ind);    
end

end