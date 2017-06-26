function A = my_Cholesky(A)
    % recursive Cholesky decomposition
    % hard to read
    [m, n] = size(A);
    chunk_size = ceil(n/2);
    if (m ~= n) error('Not a square matrix!'); end
    if (n < 20)
        A = naive_Cholesky(A);
    else
        for i = 1:chunk_size:n
            if i + chunk_size - 1 >= n
                A(i:n,i:n) = my_Cholesky(A(i:n,i:n));
            else
                A(i:i+chunk_size-1,i:i+chunk_size-1)=...
                    my_Cholesky(A(i:i+chunk_size-1,i:i+chunk_size-1));
                A(i+chunk_size:n,i:i+chunk_size-1) =...
                    (tril_solver(A(i:i+chunk_size-1,i:i+chunk_size-1),A(i+chunk_size:n,i:i+chunk_size-1)'))';
                A(i+chunk_size:n,i+chunk_size:n) = A(i+chunk_size:n,i+chunk_size:n) - A(i+chunk_size:n,i:i+chunk_size-1)*A(i+chunk_size:n,i:i+chunk_size-1)';
            end
        end
    end
end

function A = naive_Cholesky(A)
    % Cholesky decomposition based on gaxpy operation
    [m, n] = size(A);
    for j = 1 : m
        if j > 1
            A(j:n,j) = A(j:n,j) - A(j:n,1:j-1)*A(j,1:j-1)';
        end
        if (A(j,j) < 0)
            error('Not positive definite!')
        else
            A(j:n,j) = A(j:n,j)/sqrt(A(j,j));
        end
    end
end