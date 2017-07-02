clear
clc

solver_list = {'PCG','GMRES','Multi_Grid_V','Cholesky','Conjugate_Gradient','Gauss_Seidel'};

opts.k                  = 1/512;
opts.iter_num           = 512;
opts.h                  = 1/128;
opts.print_interval     = 0;
opts.Matlab_tri_solver  = 1;
opts.Matlab_Cholesky    = 0;
opts.have_some_fun      = 1;

U0 = generate_U0(128);

for solver = solver_list
    opts.subprob_solver = cell2mat(solver);
%     profile on
    [U,output] = implicit_scheme(U0,opts);
%     profsave
%     break;
end

function u0 = generate_U0 (mesh_num)
x_list = (0:mesh_num)'/mesh_num;
y_list = (0:mesh_num)'/mesh_num;
u0 = sin(pi*y_list) * sin(pi * x_list)';
end
