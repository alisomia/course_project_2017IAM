function x = CG_solver(A,b,opts)
    if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;               end
    if ~isfield(opts,'max_it');      opts.max_it  = 50;                 end
    if ~isfield(opts,'x0');          opts.x0      = zeros(size(A,2),1); end
    
    x = opts.x0;
    r = b - A * x;
    res = norm(r);
    for iter = 1 : opts.max_it
        if iter > 1
            beta = r'*r/(old_r'*old_r);
            p = r + beta * p;
        else
            p = r;
        end
        A_p = A*p;
        alpha = r'*r/(p'*A_p);
        x = x + alpha * p;
        old_r = r;
        r = r - alpha * A_p;
        if norm(A*x - b) < opts.res_tol*res
            break;
        end
    end
end