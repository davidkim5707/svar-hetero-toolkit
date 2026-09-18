function A0 = unpackA0_tvA(x, ny, lcA0, lcA0_tv, nReg)
%UNPACKA0_TVA  Reconstruct regime-stacked A0 (ny x ny x nReg) from a
%   (partial-)tvA parameter vector x with layout
%
%       x = [ a_cm (nA_cm) ; a_tv (nA_tv x nReg, regime-major) ; lambda... ]
%
%   Only the A0 block of x is read; any trailing lambda parameters are
%   ignored. This is the single unpacker shared by bvar_posterior_tvA,
%   build_theta_cache, the sampler (SECTION 12/14) and SignRestrictionCheck,
%   so every site reconstructs A0 identically.
%
%   The three behaviours collapse to one code path:
%     lcA0_tv == lcA0   -> full tvA      (a_cm empty, nA_tv = nA)
%     lcA0_tv empty/0    -> common-only   (every regime identical)
%     nReg    == 1       -> single A0     (non-tvA)
%--------------------------------------------------------------------------
    if nargin < 5 || isempty(nReg), nReg = 1; end
    if isempty(lcA0_tv), lcA0_tv = false(ny); end

    lcA0_cm = lcA0 & ~lcA0_tv;
    nA_cm   = sum(lcA0_cm(:));
    nA_tv   = sum(lcA0_tv(:));

    a_cm = x(1:nA_cm);
    A0   = zeros(ny, ny, nReg);
    for r = 1:nReg
        As = zeros(ny);
        As(lcA0_cm) = a_cm;                                  % shared across regimes
        if nA_tv > 0
            As(lcA0_tv) = x(nA_cm + (r-1)*nA_tv + (1:nA_tv)); % regime-specific block
        end
        A0(:, :, r) = As;
    end
end