% It may take several days to get the results :)
clear
clc

solution = @(x,y,t)exp(-2*pi^2*t)*sin(pi*x)*sin(pi*y);
mesh_num_list = [8,16,32,64,128,256,512];
iter_num_list = [256,512,1024,2048,4096,...
                8192,16384,65536,131072,...
                262144,524288,1048576,2097152,...
                4194304];
target_time = 1;
case_num = length(mesh_num_list);
e1 = 9999*ones(case_num,1);
e2 = 9999*ones(case_num,1);
e3 = 9999*ones(case_num,1);
e4 = 9999*ones(case_num,1);

for i = 1 : case_num
    mesh_num = mesh_num_list(i);
    acc_err= zeros(4,1);
    for iter_num = iter_num_list
        opts.iter_num           = iter_num;
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

        if opts.k/opts.h^2 < 1/8 && acc_err(1) == 0
            [U1,~] = explicit_scheme(U0,opts);
            err1 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U1) - solution(x,y,target_time))^2;
            cur_err = sqrt(int_2D(err1,[0 1],[0 1],opts))/solution_int;
            if cur_err > e1(i)
                acc_err(1) = 1;
            else
                e1(i) = cur_err;
            end
        end

        if acc_err(2) == 0
            [U2,~] = implicit_scheme(U0,opts);
            err2 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U2) - solution(x,y,target_time))^2;
            cur_err = sqrt(int_2D(err2,[0 1],[0 1],opts))/solution_int; 
            if cur_err > e2(i)
                acc_err(2) = 1;
            else
                e2(i) = cur_err;
            end
        end

        opts.theta = 1/2;
        if acc_err(3) == 0
            [U3,~] = theta_scheme(U0,opts);
            err3 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U3) - solution(x,y,target_time))^2;
            cur_err = sqrt(int_2D(err3,[0 1],[0 1],opts))/solution_int;
            if cur_err > e3(i)
                acc_err(3) = 1;
            else
                e3(i) = cur_err;
            end
        end

%         opts.theta = 1/2 - opts.h^2/(12*opts.k);
%         if opts.theta > 0 && acc_err(4) == 0
%             [U4,~] = theta_scheme(U0,opts);
%             err4 = @(x,y)(bilinear_interpolate(x,y,0,0,opts.h,U4) - solution(x,y,target_time))^2;
%             cur_err = sqrt(int_2D(err4,[0 1],[0 1],opts))/solution_int;
%             if cur_err > e4(i)
%                 acc_err(4) = 1;
%             else
%                 e4(i) = cur_err;
%             end
%         end
    end
end
loglog(1./mesh_num_list,e1,'-d',1./mesh_num_list,e2,'-v',1./mesh_num_list,e3,'-^','markerfacecolor','auto');
legend('explicit','implicit','Crank-Nicolson');
grid on

function u0 = generate_U0 (mesh_num)
x_list = (0:mesh_num)'/mesh_num;
y_list = (0:mesh_num)'/mesh_num;
u0 = sin(pi*y_list) * sin(pi * x_list)';
end
