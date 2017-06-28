function A = my_Cholesky(A,opts)
% recursive Cholesky decomposition
if ~isfield(opts,'Matlab_tri_solver'); opts.Matlab_tri_solver = 1; end
if ~isfield(opts,'threshold');   opts.threshold = 60;              end
if ~isfield(opts,'chunk_num');   opts.chunk_num =  5;              end
if ~isfield(opts,'recursive');   opts.recursive =  1;              end

[m, n] = size(A);
chunk_size = ceil(n/opts.chunk_num);
if (m ~= n) error('Not a square matrix!'); end

% solve small-scale problems directly
if (n < opts.threshold || opts.recursive == 0)
    A = naive_Cholesky(A);
else
    for i = 1:chunk_size:n
        if i + chunk_size - 1 >= n
            A(i:n,i:n) = my_Cholesky(A(i:n,i:n),opts);
        else
            cache0 = tril(my_Cholesky(A(i:i+chunk_size-1,i:i+chunk_size-1),opts));
            if opts.Matlab_tri_solver > 0
                cache1 =...
                ((cache0\A(i:i+chunk_size-1,i+chunk_size:n)))';
            else
                cache1 =...
                (tril_solver(cache0,A(i:i+chunk_size-1,i+chunk_size:n)))';
            end
            A(i:n,i:i+chunk_size-1) = [cache0;cache1];
            A(i+chunk_size:n,i+chunk_size:n) = A(i+chunk_size:n,i+chunk_size:n) - cache1*cache1';
        end
    end
end
end

function A = naive_Cholesky(A)
% Cholesky decomposition based on gaxpy operation
[m, n] = size(A);
for j = 1 : m
    if j > 1
        A(j:n,j) = A(j:n,j) - A(j:n,1:j-1)*A(j,1:j-1)' ;
    end
    if (A(j,j) < 0)
        error('Not positive definite!')
    else
        A(j:n,j) = A(j:n,j)/sqrt(A(j,j));
    end
end
end