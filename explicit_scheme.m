function [u, output] = explicit_scheme(u0,opts)
tic;

% step size in time
if ~isfield(opts,'k');                 opts.k = 1/2^10;             end
% iteration number
if ~isfield(opts,'iter_num');          opts.iter_num = 2^10;        end
% spatial step size
if ~isfield(opts,'h');                 opts.h = 1/32;             end
% draw figure every ${print_interval} iterations. (0 => off)
if ~isfield(opts,'print_interval');    opts.print_interval = 0;   end

%get mesh size
mesh_size = size(u0)-1;
x = (0:mesh_size(2))'*opts.h;
y = (0:mesh_size(1))'*opts.h;

u = u0;
progress = -1;

for iter = 1 : opts.iter_num
    u = u + opts.k * naive_laplace(u,opts.h);

    % draw figure
    if opts.print_interval > 0 && mod(iter, opts.print_interval) == 0
        surf(x,y,u);
        shading interp;
%         zlim([0 1]);
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
end

function lap = naive_laplace(u,h)
Laplace_approx = [ 0, 1, 0;...
                   1,-4, 1;...
                   0, 1, 0];
lap = conv2(u,Laplace_approx,'valid');
lap = lap / h^2;
lap = padarray(lap, [1, 1],0, 'both');
end
