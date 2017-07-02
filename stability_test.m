clear
clc

solution = @(x,y,t)exp(-2*pi^2*t)*sin(pi*x)*sin(pi*y);
mesh_num = 64;
target_time = 1;

opts.iter_num           = 512;
opts.k                  = target_time/opts.iter_num;
opts.h                  = 1/mesh_num;
opts.print_interval     = 0;
opts.theta              = 1/2;
opts.res_tol            = 1e-10;
opts.subprob_solver     = 'Conjugate_Gradient';

U0 = generate_U0(mesh_num);
x_list = (0:mesh_num)'/mesh_num;
y_list = (0:mesh_num)'/mesh_num;
solution_int = int_2D(@(x,y)solution(x,y,target_time),[0 1],[0 1],opts);


[U1,~] = explicit_scheme(U0,opts);
err1 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U1) - solution(x,y,target_time))^2;
e1 = sqrt(int_2D(err1,[0 1],[0 1],opts))/solution_int;

[U2,~] = implicit_scheme(U0,opts);
err2 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U2) - solution(x,y,target_time))^2;
e2 = sqrt(int_2D(err2,[0 1],[0 1],opts))/solution_int;

[U3,~] = theta_scheme(U0,opts);
err3 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U3) - solution(x,y,target_time))^2;
e3 = sqrt(int_2D(err3,[0 1],[0 1],opts))/solution_int;

opts.theta = 1/2 - opts.h^2/(12*opts.k);
[U4,~] = theta_scheme(U0,opts);
err4 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U4) - solution(x,y,target_time))^2;
e4 = sqrt(int_2D(err4,[0 1],[0 1],opts))/solution_int;

clc
fprintf('explicit\t:%e\n',e1);
fprintf('implicit\t:%e\n',e2);
fprintf('CN\t\t:%e\n',e3);
fprintf('surprise\t:%e\n',e4);

function u0 = generate_U0 (mesh_num)
x_list = (0:mesh_num)'/mesh_num;
y_list = (0:mesh_num)'/mesh_num;
u0 = sin(pi*y_list) * sin(pi * x_list)';
end
