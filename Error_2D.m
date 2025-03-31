  % 测试不同网格分辨率下的误差
    N_values = [50, 100, 200, 400, 800];  % 网格分辨率序列
    errors = zeros(length(N_values), 3); % 存储误差结果 [L1, L2, Linf]
    
    % 遍历所有网格分辨率
    for k = 1:length(N_values)
        % 运行仿真获取数值解
        N = N_values(k);
        WENO_2D
        
        % 计算误差
        [errors(k,1), errors(k,2), errors(k,3)] = compute_errors(u, u_exact, dx, dy);
    end
    close all
 
    orders = calc_convergence_rates(N_values, errors);
   
    function [L1, L2, Linf] = compute_errors(u, u_exact, dx, dy)
    % 计算各范数误差
    error = u - u_exact;
    
    L1 = sum(abs(error(:))) * dx * dy;
    
    L2 = sqrt(sum(error(:).^2 * dx * dy));
    
    Linf = max(abs(error(:)));

end

function orders = calc_convergence_rates(N_values, errors)
    % 计算收敛阶
    orders = zeros(length(N_values)-1, 3);
    for k = 1:length(N_values)-1
        ratio = N_values(k+1)/N_values(k);
        for m = 1:3
            orders(k,m) = log(errors(k,m)/errors(k+1,m)) / log(ratio);
        end
    end
end

 T = table();
T.N = N_values';
T.L1_Error = errors(:,1);
T.L1_Order = [nan; orders(:,1)];
T.L2_Error = errors(:,2);
T.L2_Order = [nan; orders(:,2)];
T.Linf_Error = errors(:,3);
T.Linf_Order = [nan; orders(:,3)];

disp('误差表:');
disp(T);