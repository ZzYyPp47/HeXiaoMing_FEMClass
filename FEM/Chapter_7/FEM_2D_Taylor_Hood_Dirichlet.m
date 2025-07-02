%% Local Assembly FEM for 2D example
% 计算 (u1 u2 · ∇)(u1 u2)-∇· T(u1 u2 p) = (f1 f2) on \Omega
% ∇· (u1 u2) = (0 0) on \Omega
% u1 = exp(-y) x = 0; u1 = y^2 + exp(-y) x = 1;
% u1 = x^2/16 + exp(0.25) y = -0.25; u1 = 1 y = 0;
% u2 = 2 x = 0; u2 = -2*y^3/3 + 2 x = 1;
% u2 = x/96 + 2 - pi * sin(pi * x) y = -0.25; u2 = 2 - pi * sin(pi * x) y = 0;
% f1 = -2 * mu * x^2 - 2 * mu * y ^2 - mu*exp(-y) + pi^2 * cos(pi * x) * cos(2 * pi * y);
% f2 = 4 * mu * x * y - mu * pi^3 * sin(pi * x) + 2* pi * (2 - pi * sin(pi * x)) * sin(2 * pi * y);
% mu = 1;
% real solution: [u1 u2 p] = [x^2*y^2 + exp(-y),-2*x*y^3/3 + 2 - pi * sin(pi * x),-(2 - pi * sin(pi * x)) * cos(2 * pi * y)]
% \Omega = [0,1] x [-0.25,0]
% 输入:网格单元数input_h [hx hy],是否绘图fig_flag,是否输出误差err_flag,所需误差类型err_types(用""包裹的向量)
% 输出:有限元节点处的解析解,所选err_types的误差err
function [uh,err,err_p] = FEM_2D_Taylor_Hood_Dirichlet(input_h,fig_flag,err_flag,err_types)
%% 输入网格信息
left = 0;
right = 1;
down = -0.25;
up = 0;
boundary_type = [-1 -1 -1 -1];% 边界条件 -1 -> 'D', -2 -> 'Neu', -3 -> 'Robin'

hx = input_h(1);% x 网格步长
hy = input_h(2);% y 网格步长
iteration_num = 200;% 迭代次数
tol = 10^(-8);% 最小迭代误差
dim_u = 2;% 2D
type_u = 2;% Quadratic
dim_p = 2;% 2D
type_p = 1;% Linear

[Nf_u,Nm,Nb_u,P,T,Pb_u,Tb_u,B_m_edges,B_f_edges_u,f_X_u,f_Y_u] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim_u,type_u);
[Nf_p,~,Nb_p,~,~,Pb_p,Tb_p,~,~,f_X_p,f_Y_p] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim_p,type_p);

%% 组装刚度矩阵与荷载矩阵
A1 = assembler_A(2, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, [dim_u type_u], [dim_u type_u], Tb_u, Tb_u, @func_nu, 1, 0, 1, 0);
A2 = assembler_A(2, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, [dim_u type_u], [dim_u type_u], Tb_u, Tb_u, @func_nu, 0, 1, 0, 1);
A3 = assembler_A(2, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, [dim_u type_u], [dim_u type_u], Tb_u, Tb_u, @func_nu, 1, 0, 0, 1);
A5 = assembler_A(2, Nf_p, Nf_u, Nm, P, T, Nb_p, Nb_u, [dim_p type_p], [dim_u type_u], Tb_p, Tb_u, @func_negtiveone, 0, 0, 1, 0);
A6 = assembler_A(2, Nf_p, Nf_u, Nm, P, T, Nb_p, Nb_u, [dim_p type_p], [dim_u type_u], Tb_p, Tb_u, @func_negtiveone, 0, 0, 0, 1);
b1 = assembler_b(2, Nf_u, Nm, P, T, Nb_u, [dim_u type_u], Tb_u, @func_f1, 0, 0);
b2 = assembler_b(2, Nf_u, Nm, P, T, Nb_u, [dim_u type_u], Tb_u, @func_f2, 0, 0);
O1 = sparse(Nf_p,Nf_p);
O2 = sparse(Nf_p,1);
A = [2 .* A1 + A2,A3,A5;A3',2 .* A2 + A1,A6;A5',A6',O1];
b = [b1;b2;O2];

%% 进行牛顿迭代
u0 = sparse(2 * Nf_u + Nf_p,1);
uh = Newton_iteration(iteration_num,tol,A,b,u0,B_m_edges,B_f_edges_u,Nf_u,Nf_p,Nm,P,T,Pb_u,Pb_p,Nb_u,[dim_u type_u],Tb_u,@func_boundary1,@func_boundary2,@func_real_p);
[uh_u1,uh_u2,uh_p] = deal(uh(1:Nf_u),uh(Nf_u + 1:2 * Nf_u),uh(2 * Nf_u + 1:end));

%% 求解并进行误差分析
[err,err_p] = deal([],[]);
if fig_flag == true
    show_eval(@func_real_u_1, uh_u1, Pb_u, f_X_u, f_Y_u, "Quadratic 2D FEM u1")% 调用可视化
    show_eval(@func_real_u_2, uh_u2, Pb_u, f_X_u, f_Y_u, "Quadratic 2D FEM u2")
    show_eval(@func_real_p, uh_p, Pb_p, f_X_p, f_Y_p, "Linear 2D FEM p")
end
if err_flag == true
    err_u1 = err_eval(@func_real_u_1, left, right, down, up, uh_u1, Nm, Nb_u, P, T, Pb_u, Tb_u, dim_u, type_u, input_h, err_types);
    err_u2 = err_eval(@func_real_u_2, left, right, down, up, uh_u2, Nm, Nb_u, P, T, Pb_u, Tb_u, dim_u, type_u, input_h, err_types);
    err_p = err_eval(@func_real_p, left, right, down, up, uh_p, Nm, Nb_p, P, T, Pb_p, Tb_p, dim_p, type_p, input_h, err_types);
    err = overall_err([err_u1;err_u2],err_types);
end
end