%% Local Assembly FEM for 1D example
% 计算 -(exp(x)u')' = -exp(x)(cos(x)-2sin(x)-xcos(x)-xsin(x))
% 输入:网格单元数input_h [hx hy mu],是否绘图fig_flag,是否输出误差err_flag,所需误差类型err_types(用""包裹的向量)
% 输出:有限元节点处的解析解,所选err_types的误差err
function [uh,err] = FEM_1D_Quadratic_Neumann(input_h,fig_flag,err_flag,err_types)
%% 输入网格信息
left = 0;
right = 1;
hx = input_h(1);% x 网格步长
Nx = (right - left) / hx;% x网格单元数
b_type = [-1 -2];% 边界条件 -1 -> 'D', -2 -> 'Neu', -3 -> 'Robin'
% Neu/Robin = u'+qu=p
extra_info = [-1 1];% 为 Neu/Robin 准备的额外信息
extra_info_q = [0 0];% 为 Neu/Robin 准备的额外信息(q)
extra_info_p = [0 cos(1) - sin(1)];% 为 Neu/Robin 准备的额外信息(p)


%% 输入有限元信息
Nf = 2 * Nx + 1;% x有限元节点数
Nm = Nx ;% 网格单元数
Nb = 3;% 局部基函数个数
dim = 1;% 1D
type = 2;% Quadratic

%% 输入辅助函数
real_u = @(x)x .* cos(x);
real_u_dx = @(x)cos(x) - x .* sin(x);

%% 构建P,T,Pb,Tb信息矩阵
P = left:hx:right;% 网格点坐标矩阵
T = [1:1:Nf - 1;2:1:Nf];% 网格点索引矩阵
Pb = left:hx/2:right;% 有限元点坐标矩阵
Tb = [1:2:Nf - 2;2:2:Nf - 1;3:2:Nf];% 有限元点索引矩阵
B_f_edges = [b_type;[Tb(1,1),Tb(end,end)];extra_info;extra_info_q;extra_info_p];% 有限元边界信息矩阵(第一行是边界类型,第二行是有限元边节点索引,第三/四行是额外信息(如有))

%% 组装刚度矩阵与荷载矩阵
A = assembler_A(1, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_c, 1, [], 1, []);
b = assembler_b(1, Nf, Nm, P, T, Nb, [dim type], Tb, @func_f, 0, []);

%% 使用边界条件信息矩阵调整边界条件
[A,b] = boundary_adjust(A, b, Pb, B_f_edges, @func_boundary, @func_c);% 调用调整器

%% 求解并进行误差分析
uh = A \ b;
err = [];
if fig_flag == true
    show_eval(@func_real_u, uh, Pb, "Quadratic 1D FEM")% 调用可视化
end
if err_flag == true
    err = err_eval(@func_real_u, left, right, uh, Nm, Nb, P, T, Pb, Tb, dim, type, hx, err_types);
end
end