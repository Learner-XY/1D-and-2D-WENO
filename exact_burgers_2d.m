function u_exact = exact_burgers_2d(X, Y, td ,Nx, Ny)

u0 = sin(pi*X) + sin(pi*Y);

u_exact = u0;

% 迭代参数
max_iter = 100;      % 最大迭代次数
tol = 1e-6; 

% 对每个网格点应用牛顿迭代法
for i = 1:Ny
    for j = 1:Nx
        u = u0(i, j); % 初始猜测
        for k = 1:max_iter
            % 计算隐式方程和导数
            arg_x = X(i,j) - u * td;
            arg_y = Y(i,j) - u * td;
            F = u - sin(pi*arg_x) - sin(pi*arg_y);
            dF = 1 + td*pi*(cos(pi*arg_x) + cos(pi*arg_y));

            % 牛顿迭代更新
            delta = F / dF;
            u_new = u - delta;

            % 检查收敛性
            if abs(delta) < tol
                u = u_new;
                break;
            end
            u = u_new;
        end
        u_exact(i, j) = u;
    end
end

% 绘制三维曲面图
% figure;
% surf(X, Y, u_exact);
% shading interp;
% xlabel('x');
% ylabel('y');
% zlabel('u');
% title(['2D Burgers Equation Exact Solution at t = ', num2str(t)]);
% colorbar;

% 检查收敛性
if any(isnan(u_exact(:)))
    warning('存在未收敛的点，请减小时间t或增加max_iter.');
end
end