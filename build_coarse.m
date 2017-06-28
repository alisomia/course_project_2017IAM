function [coarse_A,res_op,int_op] = build_coarse(A,usize)
% build coarse A for Multi Grid method
k = 2;
coarse_A{1} = A;
while true
    if min(usize) < 4 break; end
    res_op{k-1} = generate_I_12(usize);
    int_op{k-1} = 4 * res_op{k-1}';
    coarse_A{k} = res_op{k-1} * coarse_A{k-1} * int_op{k-1};
    usize       = (usize-1)/2;
    k           = k + 1;
end
end

function I_12 = generate_I_12(usize)
m = usize(1); n = usize(2);
kernel = [1/4, 1/2, 1/4;...
          1/2,   1, 1/2;...
          1/4, 1/2, 1/4]/4;
a = ones(m,1);
b = a;
b(1:2:m) =0;
c = a;
c(2:2:m) = 0;
d = ones(n,1);
e = d;
d(1:2:n) = 0;
e(2:2:n) = 0;
I1 = kron(spdiags(c, -1, m, m), kernel(1, 1))...
     + kron(spdiags(c, 1, m, m), kernel(1, 3))...
     + kron(spdiags(b, 0, m, m), kernel(1, 2));
I2 = kron(spdiags(c, -1, m, m), kernel(2, 1))...
     + kron(spdiags(c, 1, m, m), kernel(2, 3))...
     + kron(spdiags(b, 0, m, m), kernel(2, 2));
I3 = kron(spdiags(c, -1, m, m), kernel(3, 1))...
     + kron(spdiags(c, 1, m, m), kernel(3, 3))...
     + kron(spdiags(b, 0, m, m), kernel(3, 2));
I_12 = sparse(kron(spdiags(e, 1, n, n), I1)...
    + kron(spdiags(e, -1, n, n), I3)...
    + kron(spdiags(d, 0, n, n), I2));
index = logical(kron(d,b));
I_12 = I_12(index, :);
end