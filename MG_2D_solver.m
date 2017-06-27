function [x,res] = MG_2D_solver(A, b, mat_size, opts,coarse_A,res_op,int_op)

if ~isfield(opts,'res_tol');     opts.res_tol = 1e-6;              end
if ~isfield(opts,'max_it');      opts.max_it  = 500;               end
if ~isfield(opts,'x0');          opts.x0      = zeros(size(b));    end
if ~isfield(opts,'threshold');   opts.threshold = 8;               end
if ~isfield(opts,'Matlab_tri_solver'); opts.Matlab_tri_solver = 1; end

x = opts.x0;

% build coarse A & operators
if nargin < 5
    [coarse_A,res_op,int_op] = build_coarse(A,mat_size);
end

res0_norm = norm(A*opts.x0-b);
for iter = 1 : opts.max_it
    [x,res] = v_solver(1,b,x,mat_size,opts);
    if norm(res) < opts.res_tol*res0_norm
        break;
    end
end
    
    function [u,res] = v_solver(level,target,u0,cur_mat_size,opts)
        next_mat_size = (cur_mat_size-1)/2;
        cur_vec_size = [cur_mat_size(1)*cur_mat_size(2),1];
        next_vec_size = [next_mat_size(1)*next_mat_size(2),1];
        
        if min(cur_mat_size) < opts.threshold
            opts.x0 = u0;
            opts.max_it = 3;
            [u,res] = GS_solver(coarse_A{level},target,opts);
        else
            opts.x0 = u0;
            opts.max_it = 1;
            [u,res] = GS_solver(coarse_A{level},target,opts);
            
            coarse_res = res_op{level}*res;
%             coarse_res = restrict_operator(reshape(res,cur_mat_size));
%             coarse_res = reshape(coarse_res,next_vec_size);
            [e,~] = v_solver(level+1,coarse_res,zeros(size(coarse_res)),next_mat_size,opts);
            u = u + int_op{level}*e;
%             e = reshape(e,next_mat_size);
%             u = reshape(u,cur_mat_size) + interpolate_operator(e);
%             u = reshape(u,cur_vec_size);
            opts.max_it = 2;
            opts.x0 = u;
            [u,res] = GS_solver(coarse_A{level},target,opts);
            
        end
    end
end

function u_2h = restrict_operator(u_h)
[m,n] = size(u_h);
kernel = [1/4, 1/2, 1/4;...
          1/2,   1, 1/2;...
          1/4, 1/2, 1/4]/4;
u_h = conv2(u_h,kernel,'valid');
u_2h = u_h(1:2:m-2,1:2:n-2);
end

function u_h = interpolate_operator(u_2h)
kernel = [1/4, 1/2, 1/4;...
          1/2,   1, 1/2;...
          1/4, 1/2, 1/4];
[m,n] = size(u_2h);
u_h = zeros(2*m+1, 2*n+1);
u_h(2:2:2*m,2:2:2*n) = u_2h;
u_h = conv2(u_h,kernel,'same');
end