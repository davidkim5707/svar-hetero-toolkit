function Q = gaussian_to_Q(X)
%GAUSSIAN_TO_Q  Map n-by-n Gaussian matrix to orthogonal Q via QR.
%
%   Q = gaussian_to_Q(X)
%
%   Computes the QR decomposition of X and adjusts signs so that the
%   diagonal of R is positive.  When X has i.i.d. N(0,1) entries, the
%   resulting Q is uniformly distributed on O(n) under the Haar measure.
%
%   Reference: Rubio-Ramirez, Waggoner & Zha (2010, RES);
%              Arias, Rubio-Ramirez & Waggoner (2018, Econometrica).
%
%   INPUT
%     X : n-by-n real matrix (typically i.i.d. standard normal entries)
%
%   OUTPUT
%     Q : n-by-n orthogonal matrix  (Q'*Q = I_n)

    [Q, R] = qr(X);
    D = diag(sign(diag(R)));
    Q = Q * D;
end