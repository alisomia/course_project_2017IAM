A = full(generate_laplace_matrix(31,63,1));
coarse_A = build_coarse_A(A,[31,63]);

A = randi(3,3,5)
restrict_operator(A)
I_12 = generate_I_12([3,5]);
A0 = reshape(I_12*reshape(A,15,1),1,2)
I_21 = 4 * I_12';

interpolate_operator(A0)
reshape(I_21*reshape(A0,2,1),3,5)

function u_h = interpolate_operator(u_2h)
kernel = [1/4, 1/2, 1/4;...
          1/2,   1, 1/2;...
          1/4, 1/2, 1/4];
[m,n] = size(u_2h);
u_h = zeros(2*m+1, 2*n+1);
u_h(2:2:2*m,2:2:2*n) = u_2h;
u_h = conv2(u_h,kernel,'same');
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

function u_2h = restrict_operator(u_h)
[m,n] = size(u_h);
kernel = [1/4, 1/2, 1/4;...
          1/2,   1, 1/2;...
          1/4, 1/2, 1/4]/4;
u_h = conv2(u_h,kernel,'valid');
u_2h = u_h(1:2:m-2,1:2:n-2);
end