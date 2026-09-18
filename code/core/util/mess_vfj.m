function mess = mess_vfj(Y)
% Vats-Flegal-Jones (2019) multivariate ESS via batch means.
%   mESS = N * (|Lambda| / |Sigma_BM|)^{1/p},  batch size b = floor(N^{1/4}).
    [N, p] = size(Y);
    b = max(1, floor(N^(1/4)));
    a = floor(N / b);
    if a <= p
        warning('mess_vfj:RankShortfall', ...
            'a=%d batches <= p=%d dims; Sigma_BM rank-deficient -> NaN.', a, p);
        mess = NaN; return;
    end
    Lambda = cov(Y);
    n_use      = a * b;
    Y_use      = Y(1:n_use, :);
    batch_avg  = squeeze(mean(reshape(Y_use', p, b, a), 2))';   % a x p
    grand_avg  = mean(Y_use, 1);
    deviations = batch_avg - grand_avg;
    Sigma_BM   = (b / (a - 1)) * (deviations' * deviations);
    Lambda   = 0.5 * (Lambda   + Lambda');
    Sigma_BM = 0.5 * (Sigma_BM + Sigma_BM');
    eig_L = max(eig(Lambda),   eps);
    eig_S = max(eig(Sigma_BM), eps);
    mess = N * exp((sum(log(eig_L)) - sum(log(eig_S))) / p);
end
