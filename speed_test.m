clear
clc

addpath('./solvers');
solver_list = {'Conjugate_Gradient','Multi_Grid_V','Cholesky','Gauss_Seidel'};

opts.k                  = 1/512;
opts.iter_num           = 512;
opts.h                  = 1/64;
opts.print_interval     = 0;
opts.Matlab_tri_solver  = 1;
opts.Matlab_Cholesky    = 0;

U0 = generate_U0(64);

for solver = solver_list
    opts.subprob_solver = cell2mat(solver);
%     profile on
    [U,output] = implicit_scheme(U0,opts);
%     profsave
%     break;
end

function u0 = generate_U0 (mesh_num)
x_list = (1:mesh_num-1)'/mesh_num;
y_list = (1:mesh_num-1)'/mesh_num;
u0 = sin(pi*y_list) * sin(pi * x_list)';
end
