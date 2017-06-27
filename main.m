clear
clc

mesh_size_list = [8, 16, 32, 64, 128, 256, 512];
solver_list = {'Multi_Grid_V','PCG','Conjugate_Gradient','Cholesky','Gauss_Seidel'};

opts.k = 1/512;
opts.iter_num = 512;
opts.h = 1/128;
opts.print_interval = 0;
opts.Matlab_tri_solver = 1;
opts.Matlab_Cholesky = 1;
opts.have_some_fun = 1;

u0 = generate_u0(128);

for solver = solver_list
    opts.subprob_solver = cell2mat(solver);
    % profile on
    [u,output] = implicit_scheme(u0,opts);
    % profsave
end

function u0 = generate_u0 (mesh_num)
x = (0:mesh_num)'/mesh_num;
y = (0:mesh_num)'/mesh_num;
u0 = sin(pi*y) * sin(pi * x)';
end