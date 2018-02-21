function [U, output] = implicit_scheme(U0,opts)
tic;

% check opts. Use default value if it's not given.

% step size in time
if ~isfield(opts,'k');                  opts.k = 1/512;              end
% iteration number
if ~isfield(opts,'iter_num');           opts.iter_num = 512;         end
% spatial step size
if ~isfield(opts,'h');                  opts.h = 1/128;              end
% draw figure every ${opts.print_interval} iterations. (0 => off)
if ~isfield(opts,'print_interval');     opts.print_interval = 0;     end
% linear equation solver
if ~isfield(opts,'subprob_solver');     opts.subprob_solver = 'PCG'; end
% use built-in tril/triu equation solver ?
if ~isfield(opts,'Matlab_tri_solver');  opts.Matlab_tri_solver = 1;  end
% use built-in Cholesky factorization function ?
if ~isfield(opts,'Matlab_Cholesky');    opts.Matlab_Cholesky = 1;    end
if ~isfield(opts,'res_tol');            opts.res_tol = 1e-6;         end

mat_size = size(U0);
vec_size = [prod(mat_size),1];

% give matrix A explicitly and a handle as well
A = speye(prod(mat_size)) - opts.k * gen_Lap_2d(mat_size,opts.h);
Ax_handle = @(y)(A*y);

% Cholesky decomposition
if strcmp(opts.subprob_solver, 'Cholesky')
    if opts.Matlab_Cholesky > 0
        G = chol(A)';
    else
%         tic;
        G = tril(my_Cholesky(A,opts));
%         toc
    end
    Gt = G';
end

if strcmp(opts.subprob_solver, 'Multi_Grid_V')
%     tic;
    [coarse_A,res_op,int_op] = build_coarse(A,mat_size);
%     toc
end

if strcmp(opts.subprob_solver, 'Gauss_Seidel')
    A_struct.D_m_L = tril(A);
    A_struct.U = -triu(A,1);
end

% prepare for drawing figures
x_list = (1:mat_size(2))'*opts.h;
y_list = (1:mat_size(1))'*opts.h;


U_vec = reshape(U0,vec_size);
progress = -1;
trace.iter = zeros(opts.iter_num,1);

t1 = toc;

for iter = 1 : opts.iter_num

    switch opts.subprob_solver
        case 'Cholesky'
            if opts.Matlab_tri_solver > 0
                U_vec = G\U_vec;
                U_vec = Gt\U_vec;
            else
                opts.threshold = 50;
                opts.chunk_num = 20;
                opts.recursive = 1;
                % phase 1: G * y = b
                U_vec = tril_solver(A,U_vec,opts);
                % phase 2: G'* x = y
                U_vec = triu_solver(At,U_vec,opts);
            end
        case 'Gauss_Seidel'
            opts.x0 = U_vec;
            opts.max_it = 9999999;
            U_vec = GS_solver(A_struct,U_vec,opts);
        case 'Conjugate_Gradient'
            opts.x0 = U_vec;
            [U_vec,~,sub_output] = CG_solver(Ax_handle,U_vec,opts);
            trace.iter(iter) = sub_output.iter;
        case 'Multi_Grid_V'
            opts.x0 = U_vec;
            opts.threshold = 8;
            U_vec = MG_2D_solver(A,U_vec,mat_size,opts,coarse_A,res_op,int_op);
        otherwise
            error('Solver doesn''t exist!');
    end

    % draw figure
    if opts.print_interval > 0 && mod(iter, opts.print_interval) == 0
        mesh(x_list,y_list,reshape(U_vec,mat_size));
%         zlim([0 1]);
        drawnow;
    end

    % draw progress bar
    if (floor(100 * iter/opts.iter_num) > progress)
        progress = floor(100 * iter/opts.iter_num);
        est_time = t1+(toc-t1)*opts.iter_num/iter;
        clc
        fprintf([repmat( '=', 1, floor(progress/2)),'>', repmat('-',1,50-floor(progress/2))]);
        fprintf('\n');
        fprintf(strcat('solver \t\t:',32,opts.subprob_solver,'\n'));
        fprintf('cost time \t: %2.3f sec\nestimated time \t: %2.3f sec\n',toc,est_time);
    end
end
clc
fprintf('Done!\n');
fprintf(strcat('solver \t\t:',32,opts.subprob_solver,'\n'));
fprintf('cost time \t: %2.3f sec\n',toc);
output.cost_time = toc;
output.trace = trace;
U = reshape(U_vec,mat_size);
end
