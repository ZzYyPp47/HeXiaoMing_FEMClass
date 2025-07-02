% Nm is finite elements number
% basis_trial/test = [dim type]
% func_c need to be a function handle which can feed in (x,y).
function A = assembler_A(int_type, Nf_trial, Nf_test, Nm, P, T, Nb_trial, Nb_test, basis_trial, basis_test, Tb_trial, Tb_test, func_c, trial_dx, trial_dy, test_dx, test_dy)
    rows = cell(1, Nm); % 用于存储行索引
    cols = cell(1, Nm); % 用于存储列索引
    vals = cell(1, Nm); % 用于存储值

    for ii = 1:Nm
        S_ii = zeros(Nb_test, Nb_trial); % 初始化单元刚度矩阵
        vertices = P(:, T(1:1:size(T, 1), ii)); % 单元顶点矩阵
        for alpha = 1:Nb_trial % 计算单元刚度矩阵
            for beta = 1:Nb_test
                S_ii(beta, alpha) = get_local_integral_A(int_type, alpha, beta, basis_trial, basis_test, func_c, trial_dx, trial_dy, test_dx, test_dy, vertices);
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