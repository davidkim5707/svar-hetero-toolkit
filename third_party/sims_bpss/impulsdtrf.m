function response = impulsdtrf(By, smat, nStep, order, Q)
% impulse responses based on VAR coefficients
% This MATLAB code is translated and modified from R code originally 
% written by Karthik Sastry, part of the SVAR toolkit.
% Assumes the same model as in rfvar, except here only the By part is used.  
% smat is a square matrix of initial shock vectors.  To produce "orthogonalized
% impulse responses" it should have the property that smat'*smat=sigma, where sigma
% is the Var(u(t)) matrix and u(t) is the residual vector.  One way to get such a smat
% is to set smat=chol(sigma).  To get the smat corresponding to a different ordering,
% use smat = P' * chol(P * Sigma * P') * P, where P is a permutation matrix.  (Or equivalently,
% smat = chol(Sigma(ndx,ndx))(indx,indx), where ndx is a permutation vector and indx is its
% inverse.  E.g. ndx = [3 2 1], indx = [3,2,1], or ndx = [2, 3, 1], indx = [3,1,2].  In general
% ndx(indx) = 1:length(ndx).
% B is a neq x nvar x nlags matrix.  neq=nvar, of course, but the first index runs over 
% equations.  In response, the first index runs over variables, the second over 
% shocks (in effect, equations), the third over time.
% Code written by Christopher Sims.  This version 6/15/03.
% Inputs:
%   vout : struct, must include vout.By (lags coefficients), vout.u (residuals)
%   smat : optional, shock matrix (Cholesky or otherwise)
%   nstep: number of steps for IRF
%   order: optional, ordering vector for Cholesky decomposition
%   Q    : optional rotation matrix, default empty
%
% Output:
%   response: nvar x nshock x nstep IRF array

if nargin < 5
    Q = []; % Default: no rotation
end

if nargin < 4
    order = [];
end

if nargin < 3 || isempty(nStep)
    nStep = 40;
end

B = By;
[nvar, ~, lags] = size(B);

% Consistency check
nshock = size(smat, 2);
if size(smat,1) ~= nvar
    error('B and smat conflict on number of equations');
end

% Step 3: Initialize IR matrix
response = zeros(nvar, nshock, nStep + lags - 1);
response(:, :, lags) = smat;

% Reverse B lags
B_rev = B(:, :, lags:-1:1);
B_mat = reshape(B_rev, nvar, nvar * lags);

% Reshape for loop
response = reshape(permute(response, [1, 3, 2]), nvar * (nStep + lags - 1), nshock);

% Step 4: Recursive IR computation
ilhs = lags * nvar + (1:nvar);
irhs = 1:(nvar * lags);
for it = 1:(nStep - 1)
    response(ilhs, :) = B_mat * response(irhs, :);
    irhs = irhs + nvar;
    ilhs = ilhs + nvar;
end

% Step 5: Reshape output [nvar x horizon x nshock]
response = reshape(response, nvar, nStep + lags - 1, nshock);
response = response(:, lags:end, :);  % Discard pre-sample lags

end