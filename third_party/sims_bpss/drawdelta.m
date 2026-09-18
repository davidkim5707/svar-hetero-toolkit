function delta = drawdelta(uout, alpha, k)
% This MATLAB code is translated and modified from R code originally 
% written by Karthik Sastry, part of the SVAR toolkit.
% drawdelta: draws log scale shifts for residuals under mixture of normals
% uout: residuals (T x N)
% alpha: probability of outlier (scalar or vector)
% k: scale factor(s) for outliers (scalar or vector)

if nargin < 2
    alpha = 0.01;
end
if nargin < 3
    k = 4;
end

% Flatten residuals
u = uout(:);
T = size(uout, 1);
N = size(uout, 2);

if length(k) == 1
    % Two-component mixture
    nprob = normpdf(u);         % Standard normal density
    oprob = normpdf(u / k);     % Outlier density (scaled)
    post = (alpha * oprob) ./ (alpha * oprob + k * (1 - alpha) .* nprob);
    draw = rand(length(u), 1) < post;  % Bernoulli draw

    delta_vec = log(k) * draw;
    delta = reshape(delta_vec, T, N);

else
    % Multi-component mixture

    K = length(k);

    % Expand residuals and weights into 3D arrays
    scalemat = repmat(reshape(k, 1, 1, K), T, N, 1);      % T x N x K
    alphamat = repmat(reshape(alpha, 1, 1, K), T, N, 1);  % T x N x K
    ubig = repmat(uout, 1, 1, K);                         % T x N x K

    % Compute weighted densities
    ratio = (alphamat .* normpdf(ubig ./ scalemat)) ./ scalemat;

    % Normalize to posterior probabilities
    total = sum(ratio, 3);
    prob = ratio ./ total;

    % Cumulative probabilities for sampling
    cprob = cumsum(prob, 3);  % T x N x K

    % Random uniform draws
    rnum = rand(T, N);
    delta = zeros(T, N);

    % Log of each scale factor
    lscale = log(k);

    % Assign delta based on cumulative probability
    for iv = 1:K
        if iv == 1
            delta = delta + lscale(iv) * (rnum < cprob(:, :, iv));
        else
            delta = delta + lscale(iv) * (rnum > cprob(:, :, iv - 1) & rnum < cprob(:, :, iv));
        end
    end
end
end
