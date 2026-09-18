function delta = drawt(uout, alpha, beta)
% This MATLAB code is translated and modified from R code originally 
% written by Karthik Sastry, part of the SVAR toolkit.
% drawt: draws log standard deviation weights from inverse gamma model
% Inputs:
%   - uout: T x N matrix of residuals
%   - alpha: shape parameter (scalar or matrix compatible with uout)
%   - beta:  rate parameter (same size as alpha)
%
% Output:
%   - delta: log std devs (same size as uout)

if nargin < 2
    alpha = 4;
end
if nargin < 3
    beta = alpha;
end

% Flatten residuals
u = uout(:);  % Treat all elements symmetrically

% Prepare gamma parameters
palpha = alpha + 0.5;
pbeta = beta + 0.5 * (u .^ 2);

% Ensure palpha and pbeta are vectors of correct length
palpha = repmat(palpha, length(u), 1);
pbeta = reshape(pbeta, [], 1);

% Draw from Gamma(shape, rate)
delta_vec = gamrnd(palpha, 1 ./ pbeta);  % Gamma(shape, scale=1/rate)

% Convert to log std dev units
delta_vec = -log(delta_vec) / 2;

% Reshape to original uout shape
delta = reshape(delta_vec, size(uout));
end
