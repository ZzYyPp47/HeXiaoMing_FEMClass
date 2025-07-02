%% Local Assembly FEM for 2D example
% 计算 (u1 u2)t - ∇· T(u1 u2 p) = (f1 f2) on \Omega x T
% ∇· (u1 u2) = (0 0) on \Omega x T
% u1 = x^2*y^2 + exp(-y), t = 0;
% u2 = -2*x*y^3/3 + 2 - pi * sin(pi * x),-(2 - pi * sin(pi * x)) * cos(2 * pi * y) t = 0;
% p = -(2 - pi * sin(pi * x)) * cos(2 * pi * y) t = 0;
% f1 = -2 .* pi .* (x.^2.*y.^2+exp(-y)).*sin(2.*pi.*t)+(-2 * mu * x^2 - 2 * mu * y ^2 - mu*exp(-y) + pi^2 * cos(pi * x) * cos(2 * pi * y)).*cos(2.*pi.*t);
% f2 = -2.*pi.*(-2.*x.*y.^3./3+2-pi.*sin(pi.*x)).*sin(2.*pi.*t)+(4 * mu * x * y - mu * pi^3 * sin(pi * x) + 2* pi * (2 - pi * sin(pi * x)) * sin(2 * pi * y)).*cos(2.*pi.*t);
% mu = 1;
% real solution: [u1 u2 p] = [(x^2*y^2 + exp(-y)).*cos(2.*pi.*t),(-2*x*y^3/3 + 2 - pi * sin(pi * x)).*cos(2.*pi.*t),-(2 - pi * sin(pi * x)) * cos(2 * pi * y)*cos(2*pi*t)]
% \Omega = [0,1] x [-0.25,0] x [0,1]
% 输入:网格单元数input_h [hx hy],是否记录全部解option,是否绘图fig_flag,是否输出误差err_flag,所需误差类型err_types(用""包裹的向量)
% 输出:有限元节点处的解析解,所选err_types的误差err
function [U,err,err_p] = FEM_2D_Taylor_Hood_Dirichlet(input_h,option,fig_flag,err_flag,err_types)
%% 输入网格信息
left = 0;
right = 1;
down = -0.25;
up = 0;
t_start = 0;
t_end = 1;
boundary_type = [-1 -1 -1 -1];% 边界条件 -1 -> 'D', -2 -> 'Neu', -3 -> 'Robin'

hx = input_h(1);% x 网格步长
hy = input_h(2);% y 网格步长
ht = input_h(3);% t 网格步长
theta = input_h(4);% theta
dim_u = 2;% 2D
type_u = 2;% Quadratic
dim_p = 2;% 2D
type_p = 1;% Linear

Time = t_start:ht:t_end;% 时间坐标
[Nf_u,Nm,Nb_u,P,T,Pb_u,Tb_u,B_m_edges,B_f_edges_u,f_X_u,f_Y_u] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim_u,type_u);
[Nf_p,~,Nb_p,~,~,Pb_p,Tb_p,~,~,f_X_p,f_Y_p] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim_p,type_p);

%% 组装刚度矩阵与荷载矩阵
Me = assembler_A(2, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, [dim_u type_u], [dim_u type_u], Tb_u, Tb_u, @(x,y,t) 1, 0, 0, 0, 0, t_start);
A1 = assembler_A(2, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, [dim_u type_u], [dim_u type_u], Tb_u, Tb_u, @func_nu, 1, 0, 1, 0, t_start);
A2 = assembler_A(2, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, [dim_u type_u], [dim_u type_u], Tb_u, Tb_u, @func_nu, 0, 1, 0, 1, t_start);
A3 = assembler_A(2, Nf_u, Nf_u, Nm, P, T, Nb_u, Nb_u, [dim_u type_u], [dim_u type_u], Tb_u, Tb_u, @func_nu, 1, 0, 0, 1, t_start);
A5 = assembler_A(2, Nf_p, Nf_u, Nm, P, T, Nb_p, Nb_u, [dim_p type_p], [dim_u type_u], Tb_p, Tb_u, @(x,y,t) -1, 0, 0, 1, 0, t_start);
A6 = assembler_A(2, Nf_p, Nf_u, Nm, P, T, Nb_p, Nb_u, [dim_p type_p], [dim_u type_u], Tb_p, Tb_u, @(x,y,t) -1, 0, 0, 0, 1, t_start);
O1 = sparse(Nf_p,Nf_p);
O2 = sparse(Nf_u,Nf_p);
O3 = sparse(Nf_u,Nf_u);
O = sparse(Nf_p,1);
A = [2 .* A1 + A2,A3,A5;A3',2 .* A2 + A1,A6;A5',A6',O1];
M = [Me,O3,O2;O3,Me,O2;O2',O2',O1];
b_all = cell(1,length(Time));
parfor t_ii = 1:length(Time)
    b1 = assembler_b(2, Nf_u, Nm, P, T, Nb_u, [dim_u type_u], Tb_u, @func_f1, 0, 0, Time(t_ii));
    b2 = assembler_b(2, Nf_u, Nm, P, T, Nb_u, [dim_u type_u], Tb_u, @func_f2, 0, 0, Time(t_ii));
    b_all{t_ii} = [b1;b2;O];
end
%% 使用边界条件信息矩阵调整边界条件
u0 = [func_real_u_1(Pb_u(1,:),Pb_u(2,:),t_start,0,0)';func_real_u_2(Pb_u(1,:),Pb_u(2,:),t_start,0,0)';func_real_p(Pb_p(1,:),Pb_p(2,:),t_start,0,0)'];
[U, uh_end] = t_fdm(M, A, b_all, theta, Time, u0, option, P, T, Pb_u, Pb_p, Tb_u, Nf_u, B_f_edges_u, B_m_edges, Nb_u, @func_boundary1, @func_boundary2, @func_real_p,[],[],[],[], dim_u, type_u);

%% 求解并进行误差分析
[uh_1,uh_2,uh_p] = deal(U(1:Nf_u,:),U(Nf_u + 1:2 * Nf_u,:),U(2 * Nf_u + 1:end,:));
[err,err_p] = deal([],[]);
if fig_flag == true
    show_eval(option,@func_real_u_1,uh_1,f_X_u,f_Y_u,Time, "Quadratic 2D FEM u1")% 调用可视化
    show_eval(option,@func_real_u_2,uh_2,f_X_u,f_Y_u,Time, "Quadratic 2D FEM u2")
    show_eval(option,@func_real_p,uh_p,f_X_p,f_Y_p,Time, "Linear 2D FEM p")
end
[uh_1,uh_2,uh_p] = deal(uh_end(1:Nf_u),uh_end(Nf_u + 1:2 * Nf_u),uh_end(2 * Nf_u + 1:end));
if err_flag == true
    err_u1 = err_eval(@func_real_u_1, left, right, down, up, uh_1, Nm, Nb_u, P, T, Pb_u, Tb_u, dim_u, type_u, input_h(1:3), err_types,Time(end));
    err_u2 = err_eval(@func_real_u_2, left, right, down, up, uh_2, Nm, Nb_u, P, T, Pb_u, Tb_u, dim_u, type_u, input_h(1:3), err_types,Time(end));
    err_p = err_eval(@func_real_p, left, right, down, up, uh_p, Nm, Nb_p, P, T, Pb_p, Tb_p, dim_p, type_p, input_h(1:3), err_types,Time(end));
    err = overall_err([err_u1;err_u2],err_types);
end
end