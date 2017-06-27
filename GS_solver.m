function x = GS_solver(A,b,opts)
    if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;                  end
    if ~isfield(opts,'max_it');      opts.max_it  = 500;                   end
    if ~isfield(opts,'x0');          opts.x0      = zeros(size(A,2),1);    end
    if ~isfield(opts,'Matlab_tri_solver'); opts.Matlab_tri_solver = 1; end
    
    x = opts.x0;
    res = norm(A*x - b);
    for iter = 1 : opts.max_it
        if opts.Matlab_tri_solver > 0
            x = tril(A)\(-triu(A,1)*x+b);
        else
            x = tril_solver(A,-triu(A,1)*x+b);
        end
        if norm(A*x - b) < opts.res_tol*res
            break;
        end
    end
end