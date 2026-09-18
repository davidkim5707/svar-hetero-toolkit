function [d, fsigns, diag_all] = checkrestrictions_per_shock(restrictions, y, tol, verbose)
% Decide sign per shock independently (valid when all restrictions are per-shock).
% Inputs:
%   restrictions : cellstr (your same strings)
%   y            : (nvar × horizon × nshock) IRFs
%   tol, verbose : as before
% Outputs:
%   d       : 1 if all shocks pass under their chosen sign, else 0
%   fsigns  : 1×nshock vector of +1/-1 flips per shock
%   diag_all: struct with concatenated diagnostics

if nargin < 3 || isempty(tol), tol = 1e-12; end
if nargin < 4, verbose = false; end

[nvar, horizon, S] = size(y);
fsigns = ones(1,S);
all_pass = true;

pass_vec_cat = []; kind_cat = strings(0,1); fail_ix_cat = {};
n_sign_total = 0; n_sign_pass = 0;
n_el_total   = 0; n_el_pass   = 0;

for s = 1:S
    % ---- 1) Keep only restrictions that reference shock s
    restr_s = filter_restrictions_for_shock(restrictions, s);

    if isempty(restr_s)
        % Nothing constrains this shock; keep sign = +1
        continue
    end

    % ---- 2) Remap shock index in strings to 1, and make a (nvar×horizon×1) IRF
    restr_s1 = remap_shock_to_one(restr_s, s);
    y_s      = y(:,:,s);                    % nvar × horizon
    y_s      = reshape(y_s, [nvar, horizon, 1]);  % nvar × horizon × 1

    % ---- 3) Test original sign
    [ok_pos, diag_pos] = check_all_restrictions(y_s, restr_s1, tol, verbose);

    % ---- 4) Test flipped sign
    y_sf = -y_s;
    [ok_neg, diag_neg] = check_all_restrictions(y_sf, restr_s1, tol, verbose);

    % ---- 5) Choose
    if ok_pos || ok_neg
        if ok_pos
            chosen = diag_pos; fsigns(s) =  1;
        else
            chosen = diag_neg; fsigns(s) = -1;
        end
    else
        % Neither passes fully; pick the one with more satisfied lines (for diagnostics)
        if sum(diag_pos.pass_vec) >= sum(diag_neg.pass_vec)
            chosen = diag_pos; fsigns(s) =  1;
        else
            chosen = diag_neg; fsigns(s) = -1;
        end
        all_pass = false;
    end

    % ---- 6) Accumulate diagnostics
    pass_vec_cat = [pass_vec_cat; chosen.pass_vec(:)];
    kind_cat     = [kind_cat; chosen.kind(:)];
    fail_ix_cat  = [fail_ix_cat; chosen.fail_ix(:)];
    n_sign_total = n_sign_total + chosen.n_sign_total;
    n_sign_pass  = n_sign_pass  + chosen.n_sign_pass;
    n_el_total   = n_el_total   + chosen.n_el_total;
    n_el_pass    = n_el_pass    + chosen.n_el_pass;
end

% ---- Final outputs
d = all_pass;

diag_all.pass_vec       = pass_vec_cat;
diag_all.kind           = kind_cat;
diag_all.fail_ix        = fail_ix_cat;
diag_all.n_sign_total   = n_sign_total;
diag_all.n_sign_pass    = n_sign_pass;
diag_all.n_el_total     = n_el_total;
diag_all.n_el_pass      = n_el_pass;

end


% ===== Helpers =====

function restr_s = filter_restrictions_for_shock(restrictions, s)
% Keep only lines whose shock index equals s (works for sign & elasticity)
restr_s = {};
for k = 1:numel(restrictions)
    r = strtrim(restrictions{k});
    % Matches: y(i, t[,t2], s) OP 0
    if ~isempty(regexp(r, ['^y\(\s*\d+\s*,\s*\d+(:\s*\d+)?\s*,\s*' num2str(s) '\s*\)\s*(>=|<=|>|<)\s*0\s*$'], 'once'))
        restr_s{end+1,1} = r; %#ok<AGROW>
        continue
    end
    % Matches elasticity with same shock on both sides equal to s
    if ~isempty(regexp(r, ['^y\(\s*\d+\s*,\s*\d+(:\s*\d+)?\s*,\s*' num2str(s) '\s*\)\s*/\s*' ...
                           'y\(\s*\d+\s*,\s*\d+(:\s*\d+)?\s*,\s*' num2str(s) '\s*\)\s*(<=|>=)\s*[0-9eE+\-\.]+\s*$'], 'once'))
        restr_s{end+1,1} = r; %#ok<AGROW>
        continue
    end
end
end

function out = remap_shock_to_one(restrictions, s)
% Replace the shock index ", s)" with ", 1)" on BOTH sides (sign & elasticity).
pat_sign = [',\s*' num2str(s) '\s*\)'];
out = restrictions;
for k = 1:numel(restrictions)
    r = restrictions{k};
    % Replace shock index in any (... , s)
    r = regexprep(r, pat_sign, ', 1)');
    out{k} = r;
end
end

function [all_pass, diag] = check_all_restrictions(y, restrictions, tol, verbose)

n = numel(restrictions);
pass_vec = false(n,1);
kind     = strings(n,1);
fail_ix  = cell(n,1);

sign_needed  = false(n,1);
elastic_seen = false(n,1);
n_sign_total = 0; n_sign_pass = 0;
n_el_total   = 0; n_el_pass   = 0;

% ---- Explicit patterns (no non-capturing groups) ----
% Sign, single t
SIGN_SINGLE = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*' ...
               '(>=|<=|>|<)\s*0\s*$'];
% Sign, range t1:t2
SIGN_RANGE  = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*\)\s*' ...
               '(>=|<=|>|<)\s*0\s*$'];

% Elasticity, single t
EL_SINGLE = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*' ...
             'y\(\s*(\d+)\s*,\s*\2\s*,\s*\3\s*\)\s*' ...
             '(<=|>=)\s*' ...
             '([0-9]*\.?[0-9]+([eE][+-]?\d+)?)\s*$'];
% Elasticity, range t1:t2
EL_RANGE  = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*' ...
             'y\(\s*(\d+)\s*,\s*\2\s*:\s*\3\s*,\s*\4\s*\)\s*' ...
             '(<=|>=)\s*' ...
             '([0-9]*\.?[0-9]+([eE][+-]?\d+)?)\s*$'];

for k = 1:n
    restr = strtrim(restrictions{k});

    % ---------- SIGN: single ----------
    tok = regexp(restr, SIGN_SINGLE, 'tokens');
    if ~isempty(tok)
        kind(k) = "sign"; sign_needed(k) = true; n_sign_total = n_sign_total + 1;
        t = tok{1};
        i = str2double(t{1}); t1 = str2double(t{2}); t2 = t1; s = str2double(t{3}); op = t{4};
        vals = reshape(y(i, t1:t2, s), 1, []);
        ok = cmp_sign(vals, op, tol);
        if ok, pass_vec(k)=true; n_sign_pass = n_sign_pass + 1; else, fail_ix{k}=t1:t2; end
        continue
    end

    % ---------- SIGN: range ----------
    tok = regexp(restr, SIGN_RANGE, 'tokens');
    if ~isempty(tok)
        kind(k) = "sign"; sign_needed(k) = true; n_sign_total = n_sign_total + 1;
        t = tok{1};
        i = str2double(t{1}); t1 = str2double(t{2}); t2 = str2double(t{3}); s = str2double(t{4}); op = t{5};
        vals = reshape(y(i, t1:t2, s), 1, []);
        ok = cmp_sign(vals, op, tol);
        if ok, pass_vec(k)=true; n_sign_pass = n_sign_pass + 1; else, fail_ix{k}=t1:t2; end
        continue
    end

    % ---------- ELASTICITY: single ----------
    tok = regexp(restr, EL_SINGLE, 'tokens');
    if ~isempty(tok)
        kind(k) = "elasticity"; elastic_seen(k) = true; n_el_total = n_el_total + 1;
        t = tok{1};
        i = str2double(t{1}); tt = str2double(t{2}); s = str2double(t{3});
        j = str2double(t{4}); op = t{5}; c = str2double(t{6});
        num = reshape(y(i, tt, s), 1, []);
        den = reshape(y(j, tt, s), 1, []);
        ok = cmp_elastic(num, den, op, c, tol);
        if ok, pass_vec(k)=true; n_el_pass = n_el_pass + 1; else, fail_ix{k}=tt; end
        continue
    end

    % ---------- ELASTICITY: range ----------
    tok = regexp(restr, EL_RANGE, 'tokens');
    if ~isempty(tok)
        kind(k) = "elasticity"; elastic_seen(k) = true; n_el_total = n_el_total + 1;
        t = tok{1};
        i = str2double(t{1}); t1 = str2double(t{2}); t2 = str2double(t{3}); s = str2double(t{4});
        j = str2double(t{5}); op = t{6}; c = str2double(t{7});
        num = reshape(y(i, t1:t2, s), 1, []);
        den = reshape(y(j, t1:t2, s), 1, []);
        ok = cmp_elastic(num, den, op, c, tol);
        if ok, pass_vec(k)=true; n_el_pass = n_el_pass + 1; else, fail_ix{k}=t1:t2; end
        continue
    end

    % ---------- Unknown pattern ----------
    kind(k) = "unknown";
    pass_vec(k) = false;
    fail_ix{k} = [];
    if verbose
        fprintf('[parse] Unknown or unsupported restriction: %s\n', restr);
    end
end

% Final decision
if ~any(elastic_seen)
    all_pass = (n_sign_pass == n_sign_total) && all(pass_vec | ~sign_needed);
else
    all_pass = (n_sign_pass == n_sign_total) && (n_el_pass == n_el_total);
end

diag.pass_vec       = pass_vec;
diag.kind           = kind;
diag.fail_ix        = fail_ix;
diag.n_sign_total   = n_sign_total;
diag.n_sign_pass    = n_sign_pass;
diag.n_el_total     = n_el_total;
diag.n_el_pass      = n_el_pass;

end


% ===== helpers =====
function ok = cmp_sign(vals, op, tol)
switch op
    case '>' , ok = all(vals > 0 + tol);
    case '>=' , ok = all(vals >= 0 - tol);
    case '<' , ok = all(vals < 0 - tol);
    case '<=' , ok = all(vals <= 0 + tol);
    otherwise, ok = false;
end
end

function ok = cmp_elastic(num, den, op, c, tol)
safe = all(abs(den) > tol);
if ~safe, ok = false; return; end
ratio = num ./ den;
switch op
    case '<=' , ok = all(isfinite(ratio) & ratio <= c + tol);
    case '>=' , ok = all(isfinite(ratio) & ratio >= c - tol);
    otherwise, ok = false;
end
end
