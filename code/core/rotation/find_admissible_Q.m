function [X, Q, success, n_tried] = find_admissible_Q(check_fn, ny, max_tries, batch_size)
%FIND_ADMISSIBLE_Q  Parallel brute-force search for an admissible (X, Q).
%
%   [X, Q, success, n_tried] = find_admissible_Q(check_fn, ny, max_tries, batch_size)
%
%   Draws batches of random n-by-n Gaussian matrices in parallel (parfor),
%   maps each to an orthogonal Q via gaussian_to_Q, and returns the first
%   pair for which check_fn(Q) is true.
%
%   INPUTS
%     check_fn   : function handle  check_fn(Q) -> true/false
%     ny         : dimension n
%     max_tries  : total budget of Q candidates  (default 200,000)
%     batch_size : number of candidates per parfor batch  (default 2,000)
%
%   OUTPUTS
%     X       : n-by-n Gaussian matrix whose QR image satisfies restrictions
%     Q       : orthogonal matrix with check_fn(Q) == true
%     success : logical — true if an admissible pair was found
%     n_tried : total number of Q matrices tested

    if nargin < 3 || isempty(max_tries),  max_tries  = 200000; end
    if nargin < 4 || isempty(batch_size), batch_size = 2000;   end

    success = false;
    X = randn(ny, ny);
    Q = gaussian_to_Q(X);
    n_tried = 0;

    n_batches = ceil(max_tries / batch_size);

    for ib = 1:n_batches
        this_batch = min(batch_size, max_tries - n_tried);

        % --- Generate Q candidates in parallel ---
        pass_vec = false(this_batch, 1);
        X_batch  = cell(this_batch, 1);
        Q_batch  = cell(this_batch, 1);

        parfor ii = 1:this_batch
            Xi = randn(ny, ny);
            Qi = gaussian_to_Q(Xi);
            X_batch{ii} = Xi;
            Q_batch{ii} = Qi;
            pass_vec(ii) = check_fn(Qi);
        end

        n_tried = n_tried + this_batch;

        % --- Check if any passed ---
        idx = find(pass_vec, 1, 'first');
        if ~isempty(idx)
            X = X_batch{idx};
            Q = Q_batch{idx};
            success = true;
            fprintf('  Fallback succeeded after %d random draws.\n', n_tried);
            return
        end

        % Progress report every 10 batches
        if mod(ib, 10) == 0
            fprintf('  [Fallback] %d / %d tested (%.1f%%)\n', ...
                n_tried, max_tries, 100 * n_tried / max_tries);
        end
    end

    % warning('find_admissible_Q:notFound', ...
    %     'No admissible Q found in %d attempts.', max_tries);
end