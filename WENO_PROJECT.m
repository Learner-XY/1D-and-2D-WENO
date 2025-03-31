a = 0;
b = 4;
% N = 200;
t = 0;
td = 4;
remark = 0;

h = (b - a) / N;
x = a:h:b;  %半格点
xc = (x(1:end-1)+x(2:end)) / 2;  %整格点

epsilon = 1e-6;
CFL = 1;
U = zeros(1, N);
for i = 1:N
    U(i) = u0(xc(i));   %初值函数
end

%1是线性方程，2是burgers方程
global Choose;

Choose = 1;

if Choose == 1
    %算线性方程真解
    U_real = zeros(1, N);
    for i = 1:N
        U_real(i) = u0(xc(i) - td);
    end
else
    %用牛顿迭代算burgers方程的真解
    U_real = periodic_burgers_exact(xc, td, h);
end


%主程序
while t < td

    %下面根据方程选择最大波速
    if Choose == 1
        alpha = 1;
    else
         alpha = max(abs(U)) + 1e-10;
    end

        delta_t = CFL * h / alpha;

    if t + delta_t <= td
        t = t + delta_t;
    else
        delta_t = td - t;
        t = td;
    end
    remark = remark + 1;
    %RK3进行迭代，此时到最后一个t已经是td了
    U = TVD_RK3(U,h,delta_t);
end

%画图部分
figure
plot(xc, U, 'b*', xc, U_real,'r-');
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold')    
ylabel('u', 'FontSize', 12, 'FontWeight', 'bold')    
% title('终止时刻真解和数值解对比图', 'FontSize', 14)  
legend('数值解', '真解', 'Location', 'best')        
grid on                                              
axis tight                                           
set(gca, 'FontSize', 11)             




function v = WENO5_right(fp, i, N)
    epsilon = 1e-6;
    im2 = mod(i - 3 + N, N) + 1;  % i-2
    im1 = mod(i - 2 + N, N) + 1;  % i-1
    i0 = mod(i - 1 + N, N) + 1;   % i
    ip1 = mod(i + N, N) + 1;      % i+1
    ip2 = mod(i + 1 + N, N) + 1;  % i+2
   
    ul2 = fp(im2);
    ul1 = fp(im1);
    u = fp(i0);
    ur1 = fp(ip1);
    ur2 = fp(ip2);
    
    %计算光滑指示器beta
    beta0 = (13/12)*(u - 2*ur1 + ur2)^2 + (1/4)*(3*u - 4*ur1 + ur2)^2;
    beta1 = (13/12)*(ul1 - 2*u + ur1)^2 + (1/4)*(ul1 - ur1)^2;
    beta2 = (13/12)*(ul2 - 2*ul1 + u)^2 + (1/4)*(ul2 - 4*ul1 + 3*u)^2;

    %计算非线性权重
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

    v = w0*v0 + w1*v1 + w2*v2;
end


function v = WENO5_left(fm, i, N)
    epsilon = 1e-6;
    
    ip3 = mod(i + 2 + N, N) + 1;  % i+3
    ip2 = mod(i + 1 + N, N) + 1;  % i+2
    ip1 = mod(i + N, N) + 1;      % i+1
    i0 = mod(i - 1 + N, N) + 1;   % i
    im1 = mod(i - 2 + N, N) + 1;  % i-1

    ur3 = fm(ip3);
    ur2 = fm(ip2);
    ur1 = fm(ip1);
    u = fm(i0);
    ul1 = fm(im1);

    beta0 = (13/12) * (ur1 - 2*ur2 + ur3)^2 + (1/4) * (3*ur1 - 4*ur2 + ur3)^2;
    beta1 = (13/12) * (u - 2*ur1 + ur2)^2 + (1/4) * (u - ur2)^2;
    beta2 = (13/12) * (ul1 - 2*u + ur1)^2 + (1/4) * (ul1 - 4*u + 3*ur1)^2;


    w0 = 1/(epsilon + beta0)^2;
    w1 = 6/(epsilon + beta1)^2;
    w2 = 3/(epsilon + beta2)^2;

    sum = w0 + w1 + w2;
    w0 = w0/sum;
    w1 = w1/sum;
    w2 = w2/sum;

    v0 = (11*ur1 - 7*ur2 + 2*ur3)/6;
    v1 = (2*u + 5*ur1 - ur2)/6;
    v2 = (-ul1 + 5*u + 2*ur1)/6;

    v = w0*v0 + w1*v1 + w2*v2;
end

%做数值通量
function du = LH(u, h)
N = length(u);
global Choose;
%选择f(u),f1为线性方程,f2为burgers方程
if Choose == 1
    f = @f1;
else
    f = @f2;
end
%这里的alpha是lax-friedrichs分解中的波速，burgers方程中的
if Choose == 1
    alpha = 1;
else
    alpha = max(abs(u)) + 1e-10;
end

du = zeros(size(u));
    
    fp = 0.5*(f(u) + alpha*u);  % 正通量
    fm = 0.5*(f(u) - alpha*u);  % 负通量
    
    flux_right = zeros(1, N);

    %f(i+1/2)重构
    for i = 1:N
        % 正通量右重构，负通量左重构
        fp_plus_half = WENO5_right(fp, i, N);
        fm_plus_half = WENO5_left(fm, i, N);
        flux_right(i) = fp_plus_half + fm_plus_half;
    end
    
    for i = 1:N
        im = get_index(i, -1, N);
        du(i) = -(flux_right(i) - flux_right(im)) / h;
    end
end

function idx = get_index(i, shift, Nx)
    idx = mod(i - 1 + shift, Nx);
    if idx < 0
        idx = idx + Nx;
    end
    idx = idx + 1;
end



function U = TVD_RK3(U, h, dt)
    % 三阶TVD Runge-Kutta (SSP-RK3) 方法
    % 阶段1
    U1 = U + dt * LH(U, h);
    
    % 阶段2
    U2 = (3/4)*U + (1/4)*U1 + (1/4)*dt * LH(U1, h);
    
    % 阶段3
    U = (1/3)*U + (2/3)*U2 + (2/3)*dt * LH(U2, h);
end

%定义两个初值函数
function [b] = u0(x)
    b = (1/3)*sin(pi*x);
end

function [b] = u1(x)
    x = mod(x, 4);
    if x < 1
        b = 0;
    elseif x >= 3
        b = 0;
    else 
        b = 1;
    end
end







