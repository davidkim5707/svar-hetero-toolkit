function Sstar = compute_Sstar(y, lags)
% COMPUTE_SSTAR Compute variance-covariance matrix of univariate AR residuals
% Matches Carriero et al. (2024) approach
%
% Inputs:
%   y    - T x n matrix of data
%   lags - number of lags (p)
%
% Outputs:
%   Sstar - n x n matrix with univariate AR residual variances

[T, n] = size(y);

% Construct lagged regressors
XX = y(lags:T-1, :);
ilags = 1;
while ilags < lags
    ilags = ilags + 1;
    XX = [XX, y(lags+1-ilags:T-ilags, :)];
end

T_eff = T - lags;
X = [ones(1, T_eff); XX'];  % Include constant
Y = y(lags+1:end, :)';

% Compute univariate AR residuals for each variable
e = zeros(T_eff, n);
for i = 1:n
    beta_i = (X * X') \ (X * Y(i, :)');
    e(:, i) = Y(i, :)' - X' * beta_i;
end

% Variance-covariance matrix
Sstar = (e' * e) / T_eff;

end