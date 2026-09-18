function [d, fsigns] = checkrestrictions_per_shock_fast(SR_parsed, y, tol)
%CHECKRESTRICTIONS_PER_SHOCK_FAST  Zero-regex hot-loop sign restriction check.
%
%   [d, fsigns] = checkrestrictions_per_shock_fast(SR_parsed, y, tol)
%
%   Drop-in replacement for checkrestrictions_per_shock inside ESS hot
%   loops. Accepts pre-parsed numeric restrictions (from
%   parse_sign_restrictions) and performs pure array operations.
%
%   NO regex, NO num2str, NO int2str, NO str2double — all parsing happens
%   ONCE at setup in parse_sign_restrictions.m.
%
%   INPUTS
%     SR_parsed : output of parse_sign_restrictions(restrictions, nshock)
%     y         : (nvar × horizon × nshock) IRF tensor
%     tol       : numerical tolerance (default 1e-12)
%
%   OUTPUTS
%     d      : 1 if all shocks pass under their chosen sign flip, else 0
%     fsigns : 1×nshock vector of +1/-1 per shock
%
%   This function intentionally does NOT return the full diagnostic struct
%   of the original (pass_vec, kind, fail_ix, etc.) — those are unused
%   inside the ESS loop. For final storage, the original
%   checkrestrictions_per_shock is still called once per accepted draw.

    if nargin < 3 || isempty(tol), tol = 1e-12; end

    S        = SR_parsed.nshock;
    fsigns   = ones(1, S);
    all_pass = true;

    for s = 1:S
        has_s = SR_parsed.has_sign(s);
        has_e = SR_parsed.has_elast(s);

        if ~has_s && ~has_e
            % No restrictions on this shock — keep sign = +1
            continue
        end

        % Test original sign (+1)
        ok_pos = check_shock_block(SR_parsed.sign(s), SR_parsed.elast(s), ...
                                    has_s, has_e, y, s, +1, tol);
        % Test flipped sign (-1)
        ok_neg = check_shock_block(SR_parsed.sign(s), SR_parsed.elast(s), ...
                                    has_s, has_e, y, s, -1, tol);

        if ok_pos
            fsigns(s) = +1;
        elseif ok_neg
            fsigns(s) = -1;
        else
            fsigns(s) = +1;  % default
            all_pass  = false;
            % Early exit: as soon as one shock fails both signs, we can return.
            % (The ESS caller only cares about the binary d = all_pass.)
            d = 0;
            return
        end
    end

    d = double(all_pass);
end


function ok = check_shock_block(SS, EL, has_s, has_e, y, s, flip, tol)
%CHECK_SHOCK_BLOCK  Evaluate sign + elasticity restrictions for one shock.
%
%   Vectorized: loops over restrictions, but each restriction is handled
%   in pure numeric array ops. Early exit on first failure.

    % --- Sign restrictions ---
    if has_s
        nR = numel(SS.var);
        for k = 1:nR
            v  = SS.var(k);
            h1 = SS.h_start(k);
            h2 = SS.h_end(k);
            op = SS.op_code(k);

            % Extract and apply flip
            if h1 == h2
                vals = flip * y(v, h1, s);
            else
                vals = flip * y(v, h1:h2, s);
            end

            % Compare using op_code (integer dispatch is faster than string)
            switch op
                case  1     % >
                    if any(vals <= tol),  ok = false; return; end
                case  2     % >=
                    if any(vals <  -tol), ok = false; return; end
                case -1     % <
                    if any(vals >= -tol), ok = false; return; end
                case -2     % <=
                    if any(vals >  tol),  ok = false; return; end
                otherwise
                    ok = false; return
            end
        end
    end

    % --- Elasticity restrictions ---
    if has_e
        nR = numel(EL.num_var);
        for k = 1:nR
            nv = EL.num_var(k);
            dv = EL.den_var(k);
            h1 = EL.h_start(k);
            h2 = EL.h_end(k);
            op = EL.op_code(k);
            c  = EL.c(k);

            if h1 == h2
                num = flip * y(nv, h1, s);
                den = flip * y(dv, h1, s);
            else
                num = flip * y(nv, h1:h2, s);
                den = flip * y(dv, h1:h2, s);
            end

            % Safety: denominator must be away from zero
            if any(abs(den) <= tol), ok = false; return; end

            ratio = num ./ den;
            if ~all(isfinite(ratio)), ok = false; return; end

            switch op
                case  2     % >=
                    if any(ratio < c - tol), ok = false; return; end
                case -2     % <=
                    if any(ratio > c + tol), ok = false; return; end
                otherwise
                    ok = false; return
            end
        end
    end

    ok = true;
end