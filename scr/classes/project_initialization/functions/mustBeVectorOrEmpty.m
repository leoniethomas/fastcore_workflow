function mustBeVectorOrEmpty(x)
    if ~isvector(x) && ~isempty(x)
        error("Value must be a vector or empty.");
    end
end