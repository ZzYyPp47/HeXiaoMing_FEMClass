% Nm is finite elements number
% basis_trial/test = [dim type]
% func_c need to be a function handle which can feed in (x,y). 
function b = assembler_b(int_type, Nf_test, Nm, P, T, Nb_test, basis_test, Tb_test, func_f, test_dx, test_dy)
b = zeros(Nf_test, 1);% 初始化负载矩阵
temp = cell(1,Nm);
parfor ii = 1:Nm% 跑遍每个单元
    S = zeros(Nb_test, 1);% 初始化单元负载矩阵
    vertices = P(:,T(1:1:size(T,1),ii));% 单元顶点矩阵
    for beta = 1:1:Nb_test
        S(beta,1) = get_local_integral_b(int_type, beta, basis_test, func_f, test_dx, test_dy, vertices);
    end
    temp{ii} = S;% 将单元负载矩阵存入元胞中
end
for ii = 1:1:Nm
  b(Tb_test(:, ii),1) = b(Tb_test(:, ii),1) + temp{ii};% 将单元负载矩阵存入元胞中
end
end