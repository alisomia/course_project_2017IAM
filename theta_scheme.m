function [U, output] = theta_scheme(U0,opts)
tic;

% check opts. Use default value if it's not given.

if ~isfield(opts,'theta');             opts.theta = 1/2;           end
% step size in time
if ~isfield(opts,'k');                 opts.k = 1/128;             end
% iteration number
if ~isfield(opts,'iter_num');          opts.iter_num = 128;        end
% spatial step size
if ~isfield(opts,'h');                 opts.h = 1/32;              end
% draw figure every ${opts.print_interval} iterations. (0 => off)
if ~isfield(opts,'print_interval');    opts.print_interval = 0;    end
% linear equation solver
if ~isfield(opts,'subprob_solver');    opts.subprob_solver = 'PCG';end
if ~isfield(opts,'res_tol');           opts.res_tol = 1e-12;       end

%get mesh size
mat_size = size(U0);
vec_size = [prod(mat_size),1];

% handle for Laplace operator
Lap = gen_Lap_2d(mat_size,opts.h);

%give matrix A explicitly and a handle as well
A = speye(prod(mat_size)) - opts.theta * opts.k * Lap;
Ax_handle = @(y)(A*y);


% prepare for drawing figures
x = (0:mat_size(2))'*opts.h;
y = (0:mat_size(1))'*opts.h;


U_vec = reshape(U0,vec_size);
progress = -1;
trace.iter = zeros(opts.iter_num,1);
for iter = 1 : opts.iter_num

    switch opts.subprob_solver
        case 'Conjugate_Gradient'
            opts.x0 = U_vec;
            U_vec = CG_solver(Ax_handle,U_vec +((1-opts.theta)*opts.k)*(Lap * U_vec),opts);
        otherwise
            error('solver not exists');
    end

    % draw figure
    if opts.print_interval > 0 && mod(iter, opts.print_interval) == 0
        mesh(x,y,padarray(reshape(U_vec,m-2,n-2),[1,1],0,'both'));
        zlim([0 1]);
        drawnow;
    end

    % draw progress bar
    if (floor(100 * iter/opts.iter_num) > progress)
        progress = floor(100 * iter/opts.iter_num);
        clc
        fprintf([repmat( '=', 1, floor(progress/2)),'>', repmat('-',1,50-floor(progress/2))]);
        fprintf('\n');
    end
end

fprintf('Done!\n');
fprintf('cost time \t: %2.1f sec\n',toc);
output.cost_time = toc;
output.trace = trace;
U = reshape(U_vec,mat_size);
end
