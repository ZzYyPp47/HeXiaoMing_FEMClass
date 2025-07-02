% Nm is finite elements number
% basis_trial/test = [dim type]
% func_c need to be a function handle which can feed in (x,y).
function A = assembler_A_Newton(uh, int_type, u_dx, u_dy, Nf_trial, Nf_test, Nm, P, T, Nb_trial, Nb_test, basis_trial, basis_test, Tb_trial, Tb_test, trial_dx, trial_dy, test_dx, test_dy)
    rows = cell(1, Nm); % 用于存储行索引
    cols = cell(1, Nm); % 用于存储列索引
    vals = cell(1, Nm); % 用于存储值

    parfor ii = 1:Nm
        S_ii = zeros(Nb_test, Nb_trial); % 初始化单元刚度矩阵
        vertices = P(:, T(1:1:size(T, 1), ii)); % 单元顶点矩阵
        for alpha = 1:Nb_trial % 计算单元刚度矩阵
            for beta = 1:Nb_test
                S_ii(beta, alpha) = get_local_integral_A(uh, ii, u_dx, u_dy, P, T, Tb_trial, Nb_trial, int_type, alpha, beta, basis_trial, basis_test, trial_dx, trial_dy, test_dx, test_dy, vertices);
            end
        end
        temp_row = Tb_test(1:Nb_test, ii);
        temp_col = Tb_trial(1:Nb_trial, ii);
        [row,col] = meshgrid(temp_row,temp_col);
        [row,col] = deal(row',col');
        idx = (S_ii ~= 0);% 找到S_ii中所有非零元素及其对应索引
        rows{ii} = row(idx);
        cols{ii} = col(idx);
        vals{ii} = S_ii(idx);
    end

    % 将所有非零元素合并到全局稀疏矩阵
    rows = vertcat(rows{:});
    cols = vertcat(cols{:});
    vals = vertcat(vals{:});
    A = sparse(rows, cols, vals, Nf_test, Nf_trial);
end

%% FOR NEWTON ITERATION ONLY
% compute each integral
% basis_trial/test = [dim type]
function result = get_local_integral_A(uh, element_index, u_dx, u_dy, P, T, Tb_trial, Nb_trial, int_type, alpha, beta, basis_trial, basis_test, trial_dx, trial_dy, test_dx, test_dy, vertices)

% do 1D integral, vertices = [xmin,xmax]
if int_type == 1 
    int_func = @(x) get_element_local_basis_pvalue(uh, x, [], element_index, Nb_trial, P, T, Tb_trial, basis_trial(1), basis_trial(2), u_dx, u_dy)...
                 .* basis_local_function(x, [], basis_trial(1), basis_trial(2), alpha, trial_dx, trial_dy, vertices)...
                 .* basis_local_function(x, [], basis_test(1), basis_test(2), beta, test_dx, test_dy, vertices);
    result = gauss_1d(int_func, vertices(1), vertices(2));
end

% do 2D integral (triangle), vertices = [x1 x2 x3;y1 y2 y3]
if int_type == 2 
    int_func = @(x,y) get_element_local_basis_pvalue(uh, x, y, element_index, Nb_trial, P, T, Tb_trial, basis_trial(1), basis_trial(2), u_dx, u_dy)...
                 .* basis_local_function(x, y, basis_trial(1), basis_trial(2), alpha, trial_dx, trial_dy, vertices)...
                 .* basis_local_function(x, y, basis_test(1), basis_test(2), beta, test_dx, test_dy, vertices);
    result = gauss_2d_triangle(int_func, vertices);
end
end