function b = triu_solver(A,b,opts)
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
    b = naive_triu_solver(A,b);
else
    if opts.recursive > 0
        for i = n : -chunk_size :1
            if i - chunk_size <= 0
                b(1:i,:) = triu_solver(A(1:i,1:i),b(1:i,:),opts);
            else
                b(i-chunk_size+1:i,:) = triu_solver(A(i-chunk_size+1:i,i-chunk_size+1:i),b(i-chunk_size+1:i,:),opts);
                b(1:i-chunk_size,:) = b(1:i-chunk_size,:) - A(1:i-chunk_size,i-chunk_size+1:i)*b(i-chunk_size+1:i,:);
            end
        end
    else
       for i = n : -chunk_size :1
            if i - chunk_size <= 0
                b(1:i,:) = naive_triu_solver(A(1:i,1:i),b(1:i,:));
            else
                b(i-chunk_size+1:i,:) = naive_triu_solver(A(i-chunk_size+1:i,i-chunk_size+1:i),b(i-chunk_size+1:i,:));
                b(1:i-chunk_size,:) = b(1:i-chunk_size,:) - A(1:i-chunk_size,i-chunk_size+1:i)*b(i-chunk_size+1:i,:);
            end
       end 
    end
end
end

function u = naive_triu_solver(A_part,u)
m = size(A_part,1);
u(m,:) = u(m,:)/A_part(m,m);
for i = m-1:-1: 1
    u(i,:) = (u(i,:) - A_part(i,i+1:m)*u(i+1:m,:))/A_part(i,i);
end
end