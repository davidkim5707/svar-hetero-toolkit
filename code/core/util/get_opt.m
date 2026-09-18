function v = get_opt(opts, name, default_val)
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = default_val;
    end
end