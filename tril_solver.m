function b = tril_solver(A,b,opts)
% solves linear equation Ax = b
% A is a n*n lower triangular matrix
if nargin < 3 opts={}; end
if ~isfield(opts,'threshold');   opts.threshold = 50;               end
if ~isfield(opts,'chunk_num');   opts.chunk_num = 20;               end
if ~isfield(opts,'recursive');   opts.recursive =  1;               end


[m,n] = size(A);
chunk_size = ceil(n/opts.chunk_num);
if m~=n error('Not a square matrix!'); end

% solve small-scale problems directly
if n < opts.threshold
    b = naive_tril_solver(A,b);
else
    if opts.recursive > 0
        for i = 1 : chunk_size :n
            if i + chunk_size - 1 >= n
                b(i:n,:) = tril_solver(A(i:n,i:n),b(i:n,:),opts);
            else
                b(i:i+chunk_size-1,:) = tril_solver(A(i:i+chunk_size-1,i:i+chunk_size-1),b(i:i+chunk_size-1,:),opts);
                b(i+chunk_size:n,:) = b(i+chunk_size:n,:) - A(i+chunk_size:n,i:i+chunk_size-1)*b(i:i+chunk_size-1,:);
            end
        end
    else
        for i = 1 : chunk_size :n
            if i + chunk_size - 1 >= n
                b(i:n,:) = naive_tril_solver(A(i:n,i:n),b(i:n,:));
            else
                b(i:i+chunk_size-1,:) = naive_tril_solver(A(i:i+chunk_size-1,i:i+chunk_size-1),b(i:i+chunk_size-1,:));
                b(i+chunk_size:n,:) = b(i+chunk_size:n,:) - A(i+chunk_size:n,i:i+chunk_size-1)*b(i:i+chunk_size-1,:);
            end
        end
    end
end
end

function u = naive_tril_solver(A_part,u)
m = size(A_part,1);
u(1,:) = u(1,:)/A_part(1,1);
for i = 2 : m
    u(i,:) = (u(i,:) - A_part(i,1:i-1)*u(1:i-1,:))/A_part(i,i);
end
end
