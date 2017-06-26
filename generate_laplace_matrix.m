function lap = generate_laplace_matrix(m,n,h)
    %get Laplace approximation matrix
    Laplace_approx = [ 0, 1, 0;...
                       1,-4, 1;...
                       0, 1, 0];
    a = ones(m,1);
    Lap1 = kron(spdiags(a, -1, m, m), Laplace_approx(1, 1))...
         + kron(spdiags(a, 1, m, m), Laplace_approx(1, 3))...
         + kron(spdiags(a, 0, m, m), Laplace_approx(1, 2));
    Lap2 = kron(spdiags(a, -1, m, m), Laplace_approx(2, 1))...
         + kron(spdiags(a, 1, m, m), Laplace_approx(2, 3))...
         + kron(spdiags(a, 0, m, m), Laplace_approx(2, 2));
    Lap3 = kron(spdiags(a, -1, m, m), Laplace_approx(3, 1))...
         + kron(spdiags(a, 1, m, m), Laplace_approx(3, 3))...
         + kron(spdiags(a, 0, m, m), Laplace_approx(3, 2));
    lap = sparse(kron(spdiags(a, 1, n, n), Lap1)...
        + kron(spdiags(a, -1, n, n), Lap3)...
        + kron(spdiags(a, 0, n, n), Lap2));
    lap = lap / h^2;
end