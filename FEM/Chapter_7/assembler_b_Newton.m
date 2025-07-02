% Nm is finite elements number
% basis_trial/test = [dim type]
% func_c need to be a function handle which can feed in (x,y). 
function b = assembler_b_Newton(uh_1, uh_2, u1_dx, u1_dy, u2_dx, u2_dy, int_type, Nf_test, Nm, P, T, Tb, Nb_u, Nb_test, basis_u, basis_test, Tb_test, test_dx, test_dy)
b = zeros(Nf_test, 1);% 初始化负载矩阵
temp = cell(1,Nm);
parfor ii = 1:Nm% 跑遍每个单元
    S = zeros(Nb_test, 1);% 初始化单元负载矩阵
    vertices = P(:,T(1:1:size(T,1),ii));% 单元顶点矩阵
    for beta = 1:1:Nb_test
        S(beta,1) = get_local_integral_b(ii, uh_1, u1_dx, u1_dy, uh_2, u2_dx, u2_dy, P, T, Tb, int_type, beta, Nb_u, basis_u, basis_test, test_dx, test_dy, vertices);
    end
    temp{ii} = S;% 将单元负载矩阵存入元胞中
end
for ii = 1:1:Nm
  b(Tb_test(:, ii),1) = b(Tb_test(:, ii),1) + temp{ii};% 将单元负载矩阵存入元胞中
end
end

%% FOR NEWTON ITERATION ONLY
% compute each integral
% basis_trial/test = [dim type]
function result = get_local_integral_b(element_index, uh_1, u1_dx, u1_dy, uh_2, u2_dx, u2_dy, P, T, Tb, int_type, beta, Nb_u, basis_u, basis_test, test_dx, test_dy, vertices)

% do 1D integral, vertices = [xmin,xmax]
if int_type == 1 
    int_func = @(x) get_element_local_basis_pvalue(uh_1, x, [], element_index, Nb_u, P, T, Tb, basis_u(1), basis_u(2), u1_dx, u1_dy)...
                 .* get_element_local_basis_pvalue(uh_2, x, [], element_index, Nb_u, P, T, Tb, basis_u(1), basis_u(2), u2_dx, u2_dy)...
                 .* basis_local_function(x, [], basis_test(1), basis_test(2), beta, test_dx, test_dy, vertices);
    result = gauss_1d(int_func, vertices(1), vertices(2));
end

% do 2D integral, vertices = [x1 x2 x3;y1 y2 y3]
if int_type == 2 
    int_func = @(x,y) get_element_local_basis_pvalue(uh_1, x, y, element_index, Nb_u, P, T, Tb, basis_u(1), basis_u(2), u1_dx, u1_dy)...
                 .* get_element_local_basis_pvalue(uh_2, x, y, element_index, Nb_u, P, T, Tb, basis_u(1), basis_u(2), u2_dx, u2_dy)...
                 .* basis_local_function(x, y, basis_test(1), basis_test(2), beta, test_dx, test_dy, vertices);
    result = gauss_2d_triangle(int_func, vertices);
end
end