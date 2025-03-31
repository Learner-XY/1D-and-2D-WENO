function u = periodic_burgers_exact(x, t, h)
    % 周期为4的Burgers方程精确解
    u0 = @(x0) (1/3) * sin(pi * mod(x0, 4));  % 周期化初值
    u = zeros(size(x));

    %定义最大迭代次数和容许误差，超过此即视为存在激波
    tol = 1e-12;
    max_iter = 300;
    
    
    for i = 1:length(x)
        xi = x(i);
        x0_guess = xi; % 初始猜测

        for n = 1:length(x)
            xn = x(n);
            if abs(xn + u0(xn)*t - xi) < 2*h
                x0_guess = xn;
                break;
            elseif n == length(x)
                warning('在 x=%.3f 处初值未实现！', xi);
            end
        end


        converged = false;
        
        % 牛顿迭代求解 x0
        for iter = 1:max_iter
            F = x0_guess + u0(x0_guess)*t - xi;
            df = 1 + (pi/3) * cos(pi * mod(x0_guess,4)) * t;
            x0_new = x0_guess - F / df;
            
            if abs(x0_new - x0_guess) < tol
                converged = true;
                break;
            end
            x0_guess = x0_new;
        end
        
        if ~converged
            warning('在 x = %.3f处未收敛！', xi);
            % 激波处理：采用弱解的平均值
            u(i) = 0.5 * (u0(x0_guess + 4) + u0(x0_guess));
        else
            u(i) = u0(x0_guess);
        end
    end
end
