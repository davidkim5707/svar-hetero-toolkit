function test_log_volume_element_zero
%TEST_LOG_VOLUME_ELEMENT_ZERO  Self-test for log_volume_element_zero.m
%
%   Requires (ARRW replication package, same folder):
%       log_volume_element_zero.m
%       LogVolumeElement.m, SpheresToQ.m, SpheresRestriction.m,
%       perp.m, NumericalDerivative.m, LogAbsDet.m
%
%   Checks
%     A  single shock restricted -> function returns 0, AND the underlying ARRW
%        volume element is constant across draws (so skipping reweighting is exact).
%     B  two shocks restricted   -> q_tilde recovery reproduces Q to machine zero,
%        zero residual ~ 0, log_v finite and equal to the manual ARRW pipeline.
%     C  order dependence (sanity) -> two feasible build orders give DIFFERENT
%        per-Q volume elements (the construction measure is order dependent, which
%        is exactly why the weight must match the col_order used for sampling).

    fprintf('== test_log_volume_element_zero ==\n');
    rng(7);
    ny = 6;

    %% ---- TEST A: single-column constancy / skip ----
    ca = randn(1, ny); ca = ca / norm(ca);
    zcA.has_zero      = true;
    zcA.single_column = true;
    zcA.col_order     = [3 1 2 4 5 6];          % restricted shock (3) built first
    zcA.Z             = cell(1, ny);
    zcA.Z{3}          = ca;                       % one zero on shock 3

    v0 = log_volume_element_zero(local_build_Q(zcA, ny).', zcA, ny); % must be 0 (Q=P' sampler convention)

    nrep = 200;
    ve = zeros(nrep, 1);
    for t = 1:nrep
        Q     = local_build_Q(zcA, ny);
        ve(t) = local_direct_logve(Q, zcA, ny);   % bypass skip: raw ARRW vol element
    end
    okA = (v0 == 0) && (std(ve) < 1e-7);
    fprintf(' A (single-col): returns %.3e   raw ARRW ve mean=%.4f std=%.2e   -> %s\n', ...
            v0, mean(ve), std(ve), tf(okA));

    %% ---- TEST B: two-shock recovery / weight ----
    cb = randn(1, ny); cb = cb / norm(cb);
    cc = randn(1, ny); cc = cc / norm(cc);
    zcB.has_zero      = true;
    zcB.single_column = false;
    zcB.col_order     = [2 5 1 3 4 6];          % shocks 2 and 5 restricted, built first
    zcB.Z             = cell(1, ny);
    zcB.Z{2}          = cb;
    zcB.Z{5}          = cc;

    Q = local_build_Q(zcB, ny);
    [lv, recon_err, zero_resid] = local_check(Q, zcB, ny);
    lv2 = log_volume_element_zero(Q.', zcB, ny);  % production value (Q=P' sampler convention)
    okB = (recon_err < 1e-12) && (zero_resid < 1e-12) && isfinite(lv) && (abs(lv - lv2) < 1e-12);
    fprintf(' B (two-shock):  log_v=%.4f  |SpheresToQ(qt)-Qperm|=%.2e  zero_resid=%.2e  |lv-prod|=%.2e -> %s\n', ...
            lv, recon_err, zero_resid, abs(lv - lv2), tf(okB));

    %% ---- TEST C: order dependence (sanity) ----
    zcB2 = zcB;
    zcB2.col_order = [5 2 1 3 4 6];               % swap the two restricted shocks' order
    lvC  = log_volume_element_zero(Q.', zcB2, ny);  % Q=P' sampler convention
    okC  = abs(lv2 - lvC) > 1e-6;                  % expect a nonzero difference
    fprintf(' C (order dep.): [2 5 ..]=%.4f  [5 2 ..]=%.4f  diff=%.2e (nonzero expected) -> %s\n', ...
            lv2, lvC, abs(lv2 - lvC), tf(okC));

    if okA && okB && okC
        fprintf('ALL TESTS PASSED\n');
    else
        fprintf('SOME TESTS FAILED\n');
    end
end

% ====================== local helpers ======================

function [Zb, s, co] = local_prep(zc, ny)
    co = zc.col_order(:).';
    Zb = cell(1, ny);
    for k = 1:ny
        j = co(k);
        if j <= numel(zc.Z) && ~isempty(zc.Z{j})
            Zb{k} = zc.Z{j};
        else
            Zb{k} = zeros(0, ny);
        end
    end
    s = zeros(1, ny);
    for k = 1:ny
        s(k) = ny - (k-1) - size(Zb{k}, 1);
    end
end

function Q = local_build_Q(zc, ny)
    % reference draw with COLUMNS = shock vectors p_j (the "P" convention).
    % The sampler/production use Q = P' (rows = p_j), so callers of the
    % production function pass local_build_Q(...).' ; the direct-pipeline
    % helpers below use these columns p_j as-is.
    [Zb, s, co] = local_prep(zc, ny);
    Q  = zeros(ny, ny);
    Pb = zeros(ny, 0);
    for k = 1:ny
        Np       = perp([Pb, Zb{k}.']);
        g        = randn(s(k), 1); g = g / norm(g);
        Q(:, co(k)) = Np(:, 1:s(k)) * g;
        Pb       = [Pb, Q(:, co(k))];           %#ok<AGROW>
    end
end

function q = local_recover(Q, Zb, s, co, ny)
    q   = zeros(sum(s), 1);
    Pb  = zeros(ny, 0);
    idx = 0;
    for k = 1:ny
        Np = perp([Pb, Zb{k}.']);
        q(idx+1:idx+s(k)) = Np(:, 1:s(k)).' * Q(:, co(k));
        idx = idx + s(k);
        Pb  = [Pb, Q(:, co(k))];                %#ok<AGROW>
    end
end

function lv = local_direct_logve(Q, zc, ny)
    [Zb, s, co] = local_prep(zc, ny);
    q  = local_recover(Q, Zb, s, co, ny);
    lv = LogVolumeElement(@(x) SpheresToQ(x, Zb, ny), q, @(x) SpheresRestriction(x, Zb, ny));
end

function [lv, recon_err, zero_resid] = local_check(Q, zc, ny)
    [Zb, s, co] = local_prep(zc, ny);
    q  = local_recover(Q, Zb, s, co, ny);
    Qb    = reshape(SpheresToQ(q, Zb, ny), ny, ny);   % built columns in positions 1..ny
    Qperm = Q(:, co);                                  % Q permuted into build order
    recon_err = norm(Qb - Qperm, 'fro');
    zr = 0;
    for k = 1:ny
        if ~isempty(Zb{k})
            zr = max(zr, max(abs(Zb{k} * Q(:, co(k)))));
        end
    end
    zero_resid = zr;
    lv = LogVolumeElement(@(x) SpheresToQ(x, Zb, ny), q, @(x) SpheresRestriction(x, Zb, ny));
end

function s = tf(b)
    if b, s = 'OK'; else, s = 'FAIL'; end
end