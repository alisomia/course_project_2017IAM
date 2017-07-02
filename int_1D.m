function int = int_1D(f,x_interval,opts)
% 1D integral

if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;              end

x_start = x_interval(1); x_end = x_interval(2);
H = x_end - x_start;
I{1} = (feval(f,x_start) + feval(f,x_end))*H/2;
j = 1;
while j < 10
    sum = 0;
    for i = 1 : 2: 2^j
        sum = sum + feval(f,x_start+H*i/2^j);
    end
    I{j+1} = sum*H/2^j + I{j}/2;
    int = I{j+1};
    if abs(I{j+1}-I{j}) <= opts.res_tol*abs(I{j+1})
        break;
    end
    j = j + 1;
end
end