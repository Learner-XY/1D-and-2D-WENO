% 参数设置
% N = 200;

Nx = N;
Ny = N;
xleft = 0;
xright = 2;
yleft = 0;
yright = 2;
dx = (xright - xleft) / Nx;
dy = (yright - yleft) / Ny;
CFL = 0.5;
td = 2;
remark = 0;


% 生成网格点坐标，x和y方向均为均匀网格
[X, Y] = meshgrid(xleft + dx/2 : dx : xright, yleft + dy/2 : dy : yright);


% 初始条件
%1代表正弦函数，2代表方波函数
init = 1;

if init == 1
    u0 = @(X, Y) sin(pi * X) + sin(pi * Y);
else
    u0 = @(X, Y) double((X >= 0.5 & X <= 1.5) & (Y >= 0.5 & Y <= 1.5));
end

u = u0(X, Y);

% 通量函数 f(u) 和 g(u)
% 1代表线性方程 2代表bergers方程
Choose = 1;

if Choose  ~= 1
    f = @(u) u.^2 / 2;
    g = @(u) u.^2 / 2;
    dfdu = @(u) u;
    dgdu = @(u) u;
else
    f = @(u) u;
    g = @(u) u;
    dfdu = @(u) ones(size(u));
    dgdu = @(u) ones(size(u));
end

t = 0;

%主程序
while t < td
    % 计算最大波速 alpha，用于确定时间步长 dt
    alpha_x = max(abs(dfdu(u(:))));
    alpha_y = max(abs(dgdu(u(:))));
    alpha = max(alpha_x, alpha_y); 

    dt = CFL * min(dx, dy) / (alpha + eps);  % 时间步长，满足CFL条件

    if t + dt > td
        dt = td - t; 
    end


    [L1, ~] = LH(u, dx, dy, f, g, dfdu, dgdu);
    u1 = u + dt*L1;


    [L2, ~] = LH(u1, dx, dy, f, g, dfdu, dgdu);
    u2 = (3/4)*u + (1/4)*u1 + (1/4)*dt*L2;


    [L3, ~] = LH(u2, dx, dy, f, g, dfdu, dgdu);
    u = (1/3)*u + (2/3)*u2 + (2/3)*dt*L3;


    t = t + dt;
    remark = remark + 1;
end

%下面计算线性方程的精确解
if Choose == 1
    X_translated = mod(X - td , 2);
    Y_translated = mod(Y - td, 2);
    u_exact = u0(X_translated, Y_translated);

else 
    u_exact = exact_burgers_2d(X, Y, td, Nx, Ny);
end
    
% 结果可视化
plot_comparison(X, Y, u, u_exact, td); %真解数值解对比图

% figure('Position', [100, 100, 1200, 500]);
% subplot(1,2,2);
% contourf(X, Y, u, 20, 'LineColor', 'none');  % 绘制等高线图
% colorbar;
% xlabel('x', 'FontSize', 12, 'FontWeight', 'bold')    
% ylabel('y', 'FontSize', 12, 'FontWeight', 'bold')
% set(gca, 'FontSize', 11)             

%
% subplot(1,2,1);
% surf(X, Y, u, 'EdgeColor', 'none');
% colorbar;
% xlabel('x', 'FontSize', 12, 'FontWeight', 'bold')    
% ylabel('y', 'FontSize', 12, 'FontWeight', 'bold')
% zlabel('u', 'FontSize', 12, 'FontWeight', 'bold')
% set(gca, 'FontSize', 11)             
% view(3);


function plot_comparison(X, Y, u_num, u_exact, t)
    figure('Position', [100, 100, 1200, 500]);
    
    % 数值解
    subplot(1,2,1);
    surf(X, Y, u_num, 'EdgeColor','none');
    % title([' t=', num2str(t),'时的数值解']);
    xlabel('x', 'FontSize', 12, 'FontWeight', 'bold')    
    ylabel('y', 'FontSize', 12, 'FontWeight', 'bold')
    zlabel('u', 'FontSize', 12, 'FontWeight', 'bold')
    set(gca, 'FontSize', 11)             
    view(3); colorbar;

    z_min = min(min(u_num(:)), min(u_exact(:)));
    z_max = max(max(u_num(:)), max(u_exact(:)));

    % 设置第一个子图的 z 轴范围
    zlim([z_min, z_max]);
    
    % 精确解
    subplot(1,2,2);
    surf(X, Y, u_exact, 'EdgeColor','none');
    % title([' t=', num2str(t),'时的真解']);
    xlabel('x', 'FontSize', 12, 'FontWeight', 'bold')    
    ylabel('y', 'FontSize', 12, 'FontWeight', 'bold')
    zlabel('u', 'FontSize', 12, 'FontWeight', 'bold')
    set(gca, 'FontSize', 11)             
    view(3); colorbar;

    zlim([z_min, z_max]);
    
end


function [dudt, alpha] = LH(u, dx, dy, f, g, dfdu, dgdu)
% 计算方程的右端项，即空间导数项


[Nx, Ny] = size(u);
dudt  =  zeros(Nx, Ny); %#ok<PREALL>

% x方向导数计算
alpha_x = max(abs(dfdu(u(:))));  % x方向最大波速
dfdx = zeros(Nx, Ny);

%下面我们固定y，对x变量做一维WENO重构
%由于我们网格的产生，固定y是固定行坐标i
for i = 1:Nx
    u_slice = u(i,:);
    fp = 0.5 * (f(u_slice) + alpha_x * u_slice);
    fm = 0.5 * (f(u_slice) - alpha_x * u_slice);

    % WENO5重构正负通量
    fp_right = weno5_negative_periodic(fp);  % 正通量用左边模版重构
    fm_left =  weno5_positive_periodic(fm);  % 负通量用右边模版重构

    flux = fp_right + fm_left;

    % 计算x方向导数（f(i+1/2) - f(i-1/2))/dx
    dfdx(i,:) = (flux - circshift(flux, 1)) / dx;

end

% y方向导数计算
alpha_y = max(abs(dgdu(u(:))));
dgdy = zeros(Nx, Ny);
for j = 1:Ny
    u_slice = u(:,j)';
    gp = 0.5 * (g(u_slice) + alpha_y * u_slice);
    gm = 0.5 * (g(u_slice) - alpha_y * u_slice);

    % WENO5重构正负通量
    gp_right = weno5_negative_periodic(gp);
    gm_left = weno5_positive_periodic(gm);

    flux = gp_right + gm_left;
    dgdy(:,j) = ((flux - circshift(flux, 1)) / dy)';
end


dudt  =  -(dfdx + dgdy);
alpha = max(alpha_x, alpha_y);
end

function v_right = weno5_negative_periodic(v)

N = length(v);
v_right = zeros(size(v));
epsilon = 1e-6;

for i = 1:N
    % 周期性边界条件处理
    im2 = mod(i - 3 + N, N) + 1;  % i-2
    im1 = mod(i - 2 + N, N) + 1;  % i-1
    i0 = mod(i - 1 + N, N) + 1;   % i
    ip1 = mod(i + N, N) + 1;      % i+1
    ip2 = mod(i + 1 + N, N) + 1;  % i+2

    % 重构模板点
    ul2 = v(im2);
    ul1 = v(im1);
    u = v(i0);
    ur1 = v(ip1);
    ur2 = v(ip2);

    % 光滑指示器
    beta0 = (13/12)*(u - 2*ur1 + ur2)^2 + (1/4)*(3*u - 4*ur1 + ur2)^2;
    beta1 = (13/12)*(ul1 - 2*u + ur1)^2 + (1/4)*(ul1 - ur1)^2;
    beta2 = (13/12)*(ul2 - 2*ul1 + u)^2 + (1/4)*(ul2 - 4*ul1 + 3*u)^2;

    % 计算权重
    w0 = 3/(epsilon + beta0)^2;
    w1 = 6/(epsilon + beta1)^2;
    w2 = 1/(epsilon + beta2)^2;

    %归一化
    sum = w0 + w1 + w2;
    w0 = w0/sum;
    w1 = w1/sum;
    w2 = w2/sum;

    %重构系数
    v0 = (2*u + 5*ur1 - ur2)/6;
    v1 = (-ul1 + 5*u + 2*ur1)/6;
    v2 = (2*ul2 - 7*ul1 + 11*u)/6;


    % 加权平均
    v_right(i) = w0*v0 + w1*v1 + w2*v2;
end
end


function v_left = weno5_positive_periodic(v)


N = length(v);
v_left = zeros(size(v));
epsilon = 1e-6;

for i = 1:N

    ip3 = mod(i + 2 + N, N) + 1;  % i+3
    ip2 = mod(i + 1 + N, N) + 1;  % i+2
    ip1 = mod(i + N, N) + 1;      % i+1
    i0 = mod(i - 1 + N, N) + 1;   % i
    im1 = mod(i - 2 + N, N) + 1;  % i-1

    % 重构模板点
    ur3 = v(ip3);
    ur2 = v(ip2);
    ur1 = v(ip1);
    u = v(i0);
    ul1 = v(im1);

    %重构系数
    v0 = (11*ur1 - 7*ur2 + 2*ur3)/6;
    v1 = (2*u + 5*ur1 - ur2)/6;
    v2 = (-ul1 + 5*u + 2*ur1)/6;

    % 光滑指示器
    beta0 = (13/12) * (ur1 - 2*ur2 + ur3)^2 + (1/4) * (3*ur1 - 4*ur2 + ur3)^2;
    beta1 = (13/12) * (u - 2*ur1 + ur2)^2 + (1/4) * (u - ur2)^2;
    beta2 = (13/12) * (ul1 - 2*u + ur1)^2 + (1/4) * (ul1 - 4*u + 3*ur1)^2;

    % 计算权重
    w0 = 1/(epsilon + beta0)^2;
    w1 = 6/(epsilon + beta1)^2;
    w2 = 3/(epsilon + beta2)^2;

    sum = w0 + w1 + w2;
    w0 = w0/sum;
    w1 = w1/sum;
    w2 = w2/sum;

    % 加权平均
    v_left(i) = w0*v0 + w1*v1 + w2*v2;

end
end
