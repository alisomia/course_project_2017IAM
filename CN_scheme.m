function [u, output] = CN_scheme(u0,opts)
tic;

% check opts. Use default value if it's not given.

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

%get mesh size
mesh_size = size(u0)-1;
[m, n] = size(u0);

%give matrix A explicitly and a handle as well
A = 2 * speye((m-2)*(n-2)) - opts.k * generate_laplace_matrix(m-2,n-2,opts.h);
Laplace_approx = [ 0, 1, 0;...
                   1,-4, 1;...
                   0, 1, 0];
Ax_handle = @(y)(2 * y - opts.k * ...
    reshape(conv2(reshape(y,m-2,n-2),Laplace_approx,'same')/(opts.h^2),(m-2)*(n-2),1));

% handle for Laplace operator
lap_handle = @(y)(reshape(conv2(reshape(y,m-2,n-2),Laplace_approx,'same')/(opts.h^2),(m-2)*(n-2),1));

% prepare for drawing figures
x = (0:mesh_size(2))'*opts.h;
y = (0:mesh_size(1))'*opts.h;


u = reshape(u0(2:m-1,2:n-1),(m-2) * (n-2),1);
progress = -1;

for iter = 1 : opts.iter_num

    switch opts.subprob_solver
        case 'PCG'
            [u,~] = pcg(Ax_handle, 2 * u + opts.k * lap_handle(u),1e-10);
        otherwise
            error('solver not exists');
    end

    % draw figure
    if opts.print_interval > 0 && mod(iter, opts.print_interval) == 0
        mesh(x,y,padarray(reshape(u,m-2,n-2),[1,1],0,'both'));
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
u = padarray(reshape(u,m-2,n-2),[1,1],0,'both');
end
