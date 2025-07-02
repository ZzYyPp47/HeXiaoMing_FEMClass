%% Local Assembly FEM for 2D example
% 计算 -∇·(c∇u) = -y(1-y)(1-x-x^2/2)exp(x+y)-x(1-x/2)(-3y-y^2)exp(x+y)
% c = 1;
% u = -1.5y(1-y)exp(-1+y) x = -1;u = 0.5y(1-y)exp(1+y) x = 1;
% u = -2x(1-x/2)exp(x-1) y = -1;u = 0 y = 1;
% real solution: xy(1-x/2)(1-y)exp(x+y)
% \omega = [-1,1]^2
% 输入:网格单元数input_h [hx hy],是否绘图fig_flag,是否输出误差err_flag,所需误差类型err_types(用""包裹的向量)
% 输出:有限元节点处的解析解,所选err_types的误差err
function [uh,err] = FEM_2D_Quadratic_Dirichlet(input_h,fig_flag,err_flag,err_types)
%% 输入网格信息
left = -1;
right = 1;
down = -1;
up = 1;
boundary_type = [-1 -1 -1 -1];% 边界条件 -1 -> 'D', -2 -> 'Neu', -3 -> 'Robin'

hx = input_h(1);% x 网格步长
hy = input_h(2);% y 网格步长

dim = 2;% 2D
type = 2;% Quadratic

[Nf,Nm,Nb,P,T,Pb,Tb,B_m_edges,B_f_edges,f_X,f_Y] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim,type);

%% 组装刚度矩阵与荷载矩阵
A1 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_c, 1, 0, 1, 0);
A2 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_c, 0, 1, 0, 1);
b = assembler_b(2, Nf, Nm, P, T, Nb, [dim type], Tb, @func_f, 0, 0);
A = A1 + A2;
%% 使用边界条件信息矩阵调整边界条件
[A,b] = boundary_adjust(A,b,P,T,Pb,Tb,B_f_edges,B_m_edges,Nb,Nb,@func_boundary,@func_c,[],[],[],dim,type);% 调用调整器

%% 求解并进行误差分析
uh = A \ b;
err = [];
if fig_flag == true
    show_eval(@func_real_u, uh, Pb, f_X, f_Y, "Quadratic 2D FEM")% 调用可视化
end
if err_flag == true
    err = err_eval(@func_real_u, left, right, down, up, uh, Nm, Nb, P, T, Pb, Tb, dim, type, input_h, err_types);
end
end