function int = int_2D(f,x_interval,y_interval,opts)
% 2D integral

if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;              end

int_1D_handle = @(x)int_1D(@(y)feval(f,x,y),y_interval,opts);
int = int_1D(int_1D_handle,x_interval,opts);
end