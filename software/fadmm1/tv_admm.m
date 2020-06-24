function [x, history] = tv_admm(b, rho, lambda, x_gt)
% tv_admm  Solve total variation minimization via ADMM.
%          Use single-iteration version of Split Bregman method.
%
% [x, history] = tv_admm(b, lambda, rho)
% 
% Solves the following problem via ADMM:
% 
%   minimize  (mu / 2) * || F x - b ||_2^2 
%             + lambda * || z_x ||_1 + lambda * || z_y ||_1
%   s.t.      Fx x = z_x, Fy x = z_y
%
% where b in R^(h by w).
%
% INPUT
%  b       : Input noise-corrupted image
%  rho     : The augmented Lagrangian parameter as in Boyd, et al.
%  lambda  : The overelaxation parameter of L1 problem (always 1; use mu).
%  x_gt    : If the ground truth image is specified, compute the 
%            suboptimality of the each iterate x_t (default: [])
%
% OUTPUT
%  x       : The solution matrix that is same size as b.
%  history : A structure that contains the objective value, the primal and 
%            dual residual norms, and the tolerances for the primal and 
%            dual residual norms at each iteration.
%
% This code is modified from Boyd, et al.'s ADMM total variation demo. It
% also utilizes Split Bregman mex code by Tom Goldstein to do the G-S step.

%% Global constants and defaults
QUIET    = 1;
MAX_ITER = 1000;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;

[h, w] = size(b);

% difference operator
ex = ones(w,1);
Dx = spdiags([-ex ex], 0:1, w, w);
ey = ones(h,1);
Dy = spdiags([-ey ey], 0:1, h, h);

%% ADMM solver
x = zeros(h,w);
zx = zeros(h,w);
zy = zeros(h,w);
ux = zeros(h,w);
uy = zeros(h,w);

% the overrelaxation parameter
mu = rho*2;

if nargin < 4
    x_gt = [];
else
    obj_gt = objective(mu, b, lambda, x_gt, x_gt, x_gt);
end

if ~QUIET
    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
end

t_start = tic;

% DEBUG
debug = 0;
if debug
    hh = figure;
end

for k = 1:MAX_ITER
    % Do x-update as one sweep of G-S
    x = SplitBregmanGS(rho, mu, b, x, zx, zy, ux, uy); 
    
    % Do z-update for x and y directions
    zxold = zx;
    for i = 1 : h
        z_tmp = update_z(Dx, x(i,:)', ux(i,:)', lambda, rho);
        zx(i,:) = z_tmp';
    end
    zyold = zy;
    for j = 1 : w
        z_tmp = update_z(Dy, x(:,j), uy(:,j), lambda, rho);
        zy(:,j) = z_tmp;
    end

    % Do u-update for x and y directions
    for i = 1 : h
        u_tmp = update_u(ux(i,:)', Dx, x(i,:)', zx(i,:)');
        ux(i,:) = u_tmp';
    end
    for j = 1 : w
        u_tmp = update_u(uy(:,j), Dy, x(:,j), zy(:,j));
        uy(:,j) = u_tmp;
    end

    % diagnostics, reporting, termination checks
    val = consensus(x, Dx, zx, Dy, zy, mu, b, lambda);
    history.consensus(k) = val;
    history.objval(k)    = val.objx + val.objzx + val.objzy;
    if ~isempty(x_gt)
        history.relerr(k)    = norm(x - x_gt, 'fro') / norm(x_gt, 'fro');
        history.subopt(k)    = history.objval(k) - obj_gt;
    end

    history.r_norm(k)  = get_r_norm(Dx, Dy, x, zx, zy);
    history.s_norm(k)  = get_s_norm(rho, Dx, Dy, zx, zxold, zy, zyold);

    history.eps_pri(k) = get_eps_pri(Dx, Dy, x, zx, zy, ABSTOL, RELTOL); 
    history.eps_dual(k)= get_eps_dual(rho, Dx, Dy, ux, uy, ABSTOL, RELTOL);
    
    if debug
        figure(hh);        
        clf;
        
        subplot(121);
        imagesc(x);
        colormap(gray);
        title(sprintf('Iter: %d, Err: %f, rho = %f', ...
            k, history.relerr(k), rho));
        
        subplot(122);
        semilogy(history.objval);
        xlabel('Iterations');
        ylabel('Objective Value');
        
        pause(0.01);
    end      
    
    if ~QUIET
        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
            history.r_norm(k), history.eps_pri(k), ...
            history.s_norm(k), history.eps_dual(k), history.objval(k));
    end
    
    if (history.r_norm(k) < history.eps_pri(k) && ...
        history.s_norm(k) < history.eps_dual(k))
         break;
    end
end

if debug
    close(hh);
end

history.eTIME = toc(t_start);
history.eITER = k;

end

function val = get_r_norm(Dx, Dy, x, zx, zy)
val = 0;
[h, w] = size(x);
for i = 1 : h
    val = val + norm(Dx * x(i,:)' - zx(i,:)', 'fro');
end
for j = 1 : w
    val = val + norm(Dy * x(:,j) - zy(:,j), 'fro');
end
end

function val = get_s_norm(rho, Dx, Dy, zx, zxold, zy, zyold)
val = 0;
[h, w] = size(zx);
for i = 1 : h
    val = val + norm(-rho*Dx'*(zx(i,:)' - zxold(i,:)'), 'fro');
end
for j = 1 : w
    val = val + norm(-rho*Dy'*(zy(:,j) - zyold(:,j)), 'fro');
end
end

function val = get_eps_pri(Dx, Dy, x, zx, zy, ABSTOL, RELTOL)
val = 0;
[h, w] = size(x);
for i = 1 : h
    val = val + sqrt(w)*ABSTOL + RELTOL*max(norm(Dx*x(i,:)', 'fro'), norm(-zx(i,:)', 'fro'));
end
for j = 1 : w
    val = val + sqrt(h)*ABSTOL + RELTOL*max(norm(Dy*x(:,j), 'fro'), norm(-zy(:,j), 'fro'));
end
end

function val = get_eps_dual(rho, Dx, Dy, ux, uy, ABSTOL, RELTOL)
val = 0;
[h, w] = size(ux);
for i = 1 : h
    val = val + sqrt(w)*ABSTOL + RELTOL*norm(rho*Dx'*ux(i,:)', 'fro');
end
for j = 1 : w
    val = val + sqrt(h)*ABSTOL + RELTOL*norm(rho*Dy'*uy(:,j), 'fro');
end
end

function z = update_z(D, x, u, lambda, rho)
% z-update with relaxation
z = shrinkage(D*x + u, lambda/rho);
end

function u = update_u(u, D, x, z)
% y-update
u = u + D*x - z;
end

function obj = objective(mu, b, lambda, x, zx, zy)
% objective is Split Bregman
obj = mu*0.5*norm(x - b, 'fro')^2 + lambda*sum(sum(abs(zx))) + lambda*sum(sum(abs(zy)));
end

function y = shrinkage(a, kappa)
y = max(0, a-kappa) - max(0, -a-kappa);
end

function v = consensus(x, Dx, zx, Dy, zy, mu, b, lambda)
v.xzx = norm(Dx*x' - zx', 'fro');
v.xzy = norm(Dy*x - zy, 'fro');
v.objx = mu*0.5*norm(x - b, 'fro')^2;
v.objzx = lambda*sum(sum(abs(zx)));
v.objzy = lambda*sum(sum(abs(zy)));
end
