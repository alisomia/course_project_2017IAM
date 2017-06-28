function [x,res] = CG_solver(A,b,opts)
if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;               end
if ~isfield(opts,'max_it');      opts.max_it  = 50;                 end
if ~isfield(opts,'x0');          opts.x0      = zeros(size(b)); end

x = opts.x0;
res = b - feval(A,x);
res0_norm = norm(res);
res_norm = res0_norm;
for iter = 1 : opts.max_it
    if iter > 1
        beta = (res_norm/old_res_norm)^2;
        p = res + beta * p;
    else
        p = res;
    end
    A_p = feval(A,p);
    alpha = res_norm^2/(p'*A_p);
    x = x + alpha * p;
    old_res_norm = res_norm;
    res = res - alpha * A_p;
    res_norm = norm(res);
    if res_norm < opts.res_tol*res0_norm
        break;
    end
end
end