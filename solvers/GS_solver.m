function [x,res,output] = GS_solver(A,b,opts)
if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;              end
if ~isfield(opts,'max_it');      opts.max_it  = 100;               end
if ~isfield(opts,'x0');          opts.x0      = zeros(size(b));    end
if ~isfield(opts,'Matlab_tri_solver'); opts.Matlab_tri_solver = 1; end

x = opts.x0;
if ~isstruct(A)
    D_m_L = tril(A);
    U = -triu(A,1);
else
    D_m_L = A.D_m_L;
    U = A.U;
end

U_x = U*x;
res0_norm = norm(b-D_m_L*x+U_x);
% res_list = zeros(opts.max_it+1,1);
% res_list(1) = res0_norm;
for iter = 1 : opts.max_it
    if opts.Matlab_tri_solver > 0
        x = D_m_L\(U_x+b);
    else
        x = tril_solver(D_m_L,U_x+b);
    end
    res = -U_x;
    U_x = U*x;
    res = res+U_x;
    res_norm = norm(res);
%     res_list(iter+1) = res_norm;
    if res_norm < opts.res_tol*res0_norm
        break;
    end
end
% trace.iter_list = res_list(1:iter+1);
% output.trace = trace;
output.iter = iter;
end