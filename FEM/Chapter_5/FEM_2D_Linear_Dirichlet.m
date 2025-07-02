%% Local Assembly FEM for 2D example
% 计算 -∇· σ(u1 u2) = (f1 f2) on \Omega
% u1 = 0; u2 = 0 on \Pratial \Omega;
% f1 = -(lamda + 3mu) * (-pi^2 * sin(pi * x) * sin(pi * y)) - (lamda + mu) * ((2x - 1) * (2y - 1));
% f2 = -(lamda + 2mu) * (2 * x * (x - 1)) - (lamda + mu) * (pi^2 * cos(pi* x) * cos(pi * y)) - mu * (2 * y * (y - 1));
% lamda = 1; mu = 2;
% real solution: [u1 u2] = [sin(pi * x) * sin(pi * y),x * (x - 1) * y * (y - 1)]
% \Omega = [0,1]^2
% 输入:网格单元数input_h [hx hy],是否绘图fig_flag,是否输出误差err_flag,所需误差类型err_types(用""包裹的向量)
% 输出:有限元节点处的解析解,所选err_types的误差err
function [uh,err] = FEM_2D_Linear_Dirichlet(input_h,fig_flag,err_flag,err_types)
%% 输入网格信息
left = 0;
right = 1;
down = 0;
up = 1;
boundary_type = [-1 -1 -1 -1];% 边界条件 -1 -> 'D', -2 -> 'Neu', -3 -> 'Robin'

hx = input_h(1);% x 网格步长
hy = input_h(2);% y 网格步长

dim = 2;% 2D
type = 1;% Linear

[Nf,Nm,Nb,P,T,Pb,Tb,B_m_edges,B_f_edges,f_X,f_Y] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim,type);

%% 组装刚度矩阵与荷载矩阵
A1 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_lamda, 1, 0, 1, 0);
A2 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_nu, 1, 0, 1, 0);
A3 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_nu, 0, 1, 0, 1);
A4 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_lamda, 0, 1, 1, 0);
A5 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_nu, 1, 0, 0, 1);
A6 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_lamda, 1, 0, 0, 1);
A7 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_nu, 0, 1, 1, 0);
A8 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_lamda, 0, 1, 0, 1);
b1 = assembler_b(2, Nf, Nm, P, T, Nb, [dim type], Tb, @func_f1, 0, 0);
b2 = assembler_b(2, Nf, Nm, P, T, Nb, [dim type], Tb, @func_f2, 0, 0);
A = [A1 + 2 .* A2 + A3,A4 + A5;A6 + A7,A8 + 2 .* A3 + A2];
b = [b1;b2];
%% 使用边界条件信息矩阵调整边界条件
[A,b] = boundary_adjust(A,b,P,T,Pb,Tb,B_f_edges,B_m_edges,Nf,Nb,Nb,@func_boundary1,@func_boundary2,[],[],[],[],dim,type);% 调用调整器

%% 求解并进行误差分析
uh = A \ b;
[uh_1,uh_2] = deal(uh(1:1:Nf),uh(Nf + 1:1:end));
err = [];
if fig_flag == true
    show_eval(@func_real_u_1, uh_1, Pb, f_X, f_Y, "Linear 2D FEM u1")% 调用可视化
    show_eval(@func_real_u_2, uh_2, Pb, f_X, f_Y, "Linear 2D FEM u2")
end
if err_flag == true
    err_u1 = err_eval(@func_real_u_1, left, right, down, up, uh_1, Nm, Nb, P, T, Pb, Tb, dim, type, input_h, err_types);
    err_u2 = err_eval(@func_real_u_2, left, right, down, up, uh_2, Nm, Nb, P, T, Pb, Tb, dim, type, input_h, err_types);
    err = overall_err([err_u1;err_u2],err_types);
end
end