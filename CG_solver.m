function x = CG_solver(A,b,opts)
    if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;               end
    if ~isfield(opts,'max_it');      opts.max_it  = 50;                 end
    if ~isfield(opts,'x0');          opts.x0      = zeros(size(A,2),1); end
    
    x = opts.x0;
    r = b - feval(A,x);
    res0 = norm(r);
    res = res0;
    for iter = 1 : opts.max_it
        if iter > 1
            beta = (res/old_res)^2;
            p = r + beta * p;
        else
            p = r;
        end
        A_p = feval(A,p);
        alpha = res^2/(p'*A_p);
        x = x + alpha * p;
        old_res = res;
        r = r - alpha * A_p;
        res = norm(r);
        if res < opts.res_tol*res0
            break;
        end
    end
end