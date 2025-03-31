%计算误差表
% 定义不同网格分辨率
N_list = [50, 100, 200, 400, 800]; % 示例网格数
num_cases = length(N_list);

% 预存储误差
error_L1 = zeros(num_cases, 1);
error_L2 = zeros(num_cases, 1);
error_Linf = zeros(num_cases, 1);

for k = 1:num_cases
    N = N_list(k);
    WENO_PROJECT

    u_exact = U_real;
    u_num = U; 
    
    % 计算误差
    error = u_num - u_exact;
    
    % L1误差
    error_L1(k) = sum(abs(error)) * h;
    
    % L2误差
    error_L2(k) = sqrt(sum(error.^2) * h);
    
    % L∞误差
    error_Linf(k) = max(abs(error));
end

% 计算收敛阶
order_L1 = zeros(num_cases-1, 1);
order_L2 = zeros(num_cases-1, 1);
order_Linf = zeros(num_cases-1, 1);

for k = 1:num_cases-1
    ratio = N_list(k+1)/N_list(k); % 网格加密比
    order_L1(k) = log(error_L1(k)/error_L1(k+1)) / log(ratio);
    order_L2(k) = log(error_L2(k)/error_L2(k+1)) / log(ratio);
    order_Linf(k) = log(error_Linf(k)/error_Linf(k+1)) / log(ratio);
end

% 构建并显示表格
T = table();
T.N = N_list';
T.L1_Error = error_L1;
T.L1_Order = [nan; order_L1];
T.L2_Error = error_L2;
T.L2_Order = [nan; order_L2];
T.Linf_Error = error_Linf;
T.Linf_Order = [nan; order_Linf];

disp('误差表:');
disp(T);

