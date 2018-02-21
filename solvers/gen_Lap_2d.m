function Lap = gen_Lap_2d(mat_size,h)
%get 2-D Laplace approximation matrix
m = mat_size(1);
n = mat_size(2);
Lap = kron(gen_Lap_1d(n),speye(m))+kron(speye(n),gen_Lap_1d(m));
Lap = Lap/h^2;
    function Lap_ = gen_Lap_1d(n)
        Lap_ = spdiags(-2*ones(n,1),0,n,n)...
            + spdiags(ones(n,1),1,n,n)...
            + spdiags(ones(n,1),-1,n,n);
    end
end