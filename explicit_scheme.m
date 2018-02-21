function [U, output] = explicit_scheme(U0,opts)
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
mat_size = size(U0);
vec_size = [prod(mat_size),1];

x = (1:mat_size(2))'*opts.h;
y = (1:mat_size(1))'*opts.h;

U_vec = reshape(U0,vec_size);
progress = -1;
Lap = gen_Lap_2d(mat_size,opts.h);
for iter = 1 : opts.iter_num
    U_vec = U_vec + opts.k*(Lap*U_vec);

    % draw figure
    if opts.print_interval > 0 && mod(iter, opts.print_interval) == 0
        surf(x,y,reshape(U_vec,mat_size));
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
U = reshape(U_vec,mat_size);
end