function u = exact_burgers(x_grid, t)

    
    %最大迭代次数与容许误差
    max_iter = 100;
    tol = 1e-12;
    
    u = zeros(size(x_grid)); %初始化解数组

    for i = 1:length(x_grid)
        xi = x_grid(i);
        u_guess = (1/3)*sin(pi * xi); %初始猜测: t=0时的值
        
        %牛顿迭代
        for iter = 1:max_iter
            F = u_guess - (1/3)*sin(pi*(xi - u_guess * t));
            dF = 1 + (1/3)*pi * t * cos(pi * (xi - u_guess * t));
            delta = F / dF;
            u_guess = u_guess - delta;
            
            if abs(delta) < tol
                break;
            end
        end
        
        if iter == max_iter
            warning('在 x=%.3f 处未收敛！', xi);
        end
        
        u(i) = u_guess;
    end
end
