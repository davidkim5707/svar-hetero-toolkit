function SR = parse_sign_restrictions(restrictions, nshock)
%PARSE_SIGN_RESTRICTIONS  Parse sign restriction strings into numeric arrays.
%
%   SR = parse_sign_restrictions(restrictions, nshock)
%
%   Converts a cell array of restriction strings (e.g. 'y(2,1:6,1) < 0',
%   'y(1,1,2)/y(3,1,2) <= 0.0258') into pre-computed numeric arrays
%   partitioned by shock. This is called ONCE at setup and the result is
%   then passed to the hot-loop checker (checkrestrictions_per_shock_fast).
%
%   This eliminates ALL regex / num2str / str2double calls from the hot
%   loop, which were responsible for ~40% of runtime in profiling.
%
%   INPUT
%     restrictions : cell array of strings
%     nshock       : number of shocks (typically ny)
%
%   OUTPUT
%     SR : struct with fields indexed by shock s = 1..nshock
%       SR.sign(s)   : struct-of-arrays of sign restrictions for shock s
%       SR.elast(s)  : struct-of-arrays of elasticity restrictions for shock s
%       SR.has_sign(s), SR.has_elast(s)   : logicals
%       SR.n_sign(s),  SR.n_elast(s)      : counts
%
%   OP CODE CONVENTION (integer, faster than string compare)
%     For sign restrictions (y op 0):
%       +2 = '>='    +1 = '>'     -1 = '<'     -2 = '<='
%     For elasticity restrictions (num/den op c):
%       +2 = '>='    -2 = '<='

    if nargin < 2, nshock = 10; end

    SR.nshock    = nshock;
    SR.has_sign  = false(nshock, 1);
    SR.has_elast = false(nshock, 1);
    SR.n_sign    = zeros(nshock, 1);
    SR.n_elast   = zeros(nshock, 1);

    % --- Regex patterns ---
    SIGN_SINGLE = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*' ...
                   '(>=|<=|>|<)\s*0\s*$'];
    SIGN_RANGE  = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*\)\s*' ...
                   '(>=|<=|>|<)\s*0\s*$'];
    EL_SINGLE   = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*' ...
                   'y\(\s*(\d+)\s*,\s*\2\s*,\s*\3\s*\)\s*' ...
                   '(<=|>=)\s*([0-9]*\.?[0-9]+([eE][+-]?\d+)?)\s*$'];
    EL_RANGE    = ['^y\(\s*(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*\)\s*/\s*' ...
                   'y\(\s*(\d+)\s*,\s*\2\s*:\s*\3\s*,\s*\4\s*\)\s*' ...
                   '(<=|>=)\s*([0-9]*\.?[0-9]+([eE][+-]?\d+)?)\s*$'];

    % ------------------------------------------------------------------
    % Initialize per-shock accumulators as column vectors (FIX: no struct
    % array growth — we append to numeric columns directly). This avoids
    % MATLAB's "Subscripted assignment between dissimilar structures" when
    % two entry types (sign vs elasticity) are mixed.
    % ------------------------------------------------------------------
    sign_var   = cell(nshock, 1);
    sign_h1    = cell(nshock, 1);
    sign_h2    = cell(nshock, 1);
    sign_op    = cell(nshock, 1);
    sign_opstr = cell(nshock, 1);

    elast_nv   = cell(nshock, 1);
    elast_dv   = cell(nshock, 1);
    elast_h1   = cell(nshock, 1);
    elast_h2   = cell(nshock, 1);
    elast_op   = cell(nshock, 1);
    elast_c    = cell(nshock, 1);

    for s = 1:nshock
        sign_var{s}   = [];  sign_h1{s}  = [];  sign_h2{s}  = [];
        sign_op{s}    = [];  sign_opstr{s} = {};
        elast_nv{s}   = [];  elast_dv{s} = [];
        elast_h1{s}   = [];  elast_h2{s} = [];  elast_op{s} = [];  elast_c{s}  = [];
    end

    if isempty(restrictions)
        SR.sign  = build_empty_sign(nshock);
        SR.elast = build_empty_elast(nshock);
        return
    end

    % --- Parse each restriction into numeric arrays ---
    for k = 1:numel(restrictions)
        r = strtrim(restrictions{k});
        matched = false;

        % SIGN_SINGLE
        tok = regexp(r, SIGN_SINGLE, 'tokens', 'once');
        if ~isempty(tok)
            var_i  = str2double(tok{1});
            h1     = str2double(tok{2});
            s      = str2double(tok{3});
            op_str = tok{4};
            sign_var{s}(end+1,1)   = var_i;
            sign_h1{s}(end+1,1)    = h1;
            sign_h2{s}(end+1,1)    = h1;
            sign_op{s}(end+1,1)    = op_to_code(op_str);
            sign_opstr{s}{end+1,1} = op_str;
            matched = true;
        end

        % SIGN_RANGE
        if ~matched
            tok = regexp(r, SIGN_RANGE, 'tokens', 'once');
            if ~isempty(tok)
                var_i  = str2double(tok{1});
                h1     = str2double(tok{2});
                h2     = str2double(tok{3});
                s      = str2double(tok{4});
                op_str = tok{5};
                sign_var{s}(end+1,1)   = var_i;
                sign_h1{s}(end+1,1)    = h1;
                sign_h2{s}(end+1,1)    = h2;
                sign_op{s}(end+1,1)    = op_to_code(op_str);
                sign_opstr{s}{end+1,1} = op_str;
                matched = true;
            end
        end

        % EL_SINGLE
        if ~matched
            tok = regexp(r, EL_SINGLE, 'tokens', 'once');
            if ~isempty(tok)
                num_v  = str2double(tok{1});
                h1     = str2double(tok{2});
                s      = str2double(tok{3});
                den_v  = str2double(tok{4});
                op_str = tok{5};
                c      = str2double(tok{6});
                elast_nv{s}(end+1,1) = num_v;
                elast_dv{s}(end+1,1) = den_v;
                elast_h1{s}(end+1,1) = h1;
                elast_h2{s}(end+1,1) = h1;
                elast_op{s}(end+1,1) = op_to_code(op_str);
                elast_c{s}(end+1,1)  = c;
                matched = true;
            end
        end

        % EL_RANGE
        if ~matched
            tok = regexp(r, EL_RANGE, 'tokens', 'once');
            if ~isempty(tok)
                num_v  = str2double(tok{1});
                h1     = str2double(tok{2});
                h2     = str2double(tok{3});
                s      = str2double(tok{4});
                den_v  = str2double(tok{5});
                op_str = tok{6};
                c      = str2double(tok{7});
                elast_nv{s}(end+1,1) = num_v;
                elast_dv{s}(end+1,1) = den_v;
                elast_h1{s}(end+1,1) = h1;
                elast_h2{s}(end+1,1) = h2;
                elast_op{s}(end+1,1) = op_to_code(op_str);
                elast_c{s}(end+1,1)  = c;
                matched = true;
            end
        end

        if ~matched
            warning('parse_sign_restrictions:unknown', ...
                    'Unknown restriction pattern: %s', r);
        end
    end

    % --- Finalize ---
    SR.sign  = build_empty_sign(nshock);
    SR.elast = build_empty_elast(nshock);
    for s = 1:nshock
        SR.sign(s).var     = sign_var{s};
        SR.sign(s).h_start = sign_h1{s};
        SR.sign(s).h_end   = sign_h2{s};
        SR.sign(s).op_code = sign_op{s};
        SR.sign(s).op_str  = sign_opstr{s};

        SR.elast(s).num_var = elast_nv{s};
        SR.elast(s).den_var = elast_dv{s};
        SR.elast(s).h_start = elast_h1{s};
        SR.elast(s).h_end   = elast_h2{s};
        SR.elast(s).op_code = elast_op{s};
        SR.elast(s).c       = elast_c{s};

        SR.n_sign(s)  = numel(sign_var{s});
        SR.n_elast(s) = numel(elast_nv{s});
        SR.has_sign(s)  = SR.n_sign(s)  > 0;
        SR.has_elast(s) = SR.n_elast(s) > 0;
    end
end


function code = op_to_code(op_str)
    switch op_str
        case '>=', code =  2;
        case '>',  code =  1;
        case '<',  code = -1;
        case '<=', code = -2;
        otherwise, code =  0;
    end
end


function out = build_empty_sign(n)
    template = struct('var',[],'h_start',[],'h_end',[],'op_code',[],'op_str',{{}});
    out = repmat(template, n, 1);
end


function out = build_empty_elast(n)
    template = struct('num_var',[],'den_var',[],'h_start',[],'h_end',[], ...
                      'op_code',[],'c',[]);
    out = repmat(template, n, 1);
end