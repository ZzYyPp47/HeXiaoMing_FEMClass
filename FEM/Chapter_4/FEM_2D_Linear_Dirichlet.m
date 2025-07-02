%% Local Assembly FEM for 2D example
% 计算 ut-∇·(c∇u) = -3exp(x+y+t)
% c = 2;
% u = exp(x+y) t = 0 & in \Omega;u = exp(y+t) x = 0;
% u = exp(2+y+t) x = 2;u = exp(x+t) y = 0;u = exp(x+1+t) y = 1;
% real solution: exp(x+y+t)
% \omega x T = [0,2] x [0,1] x [0,1]
% 注意:此程序会记录整个\Omega x T 上的有限元节点处的解,而不是仅终了时刻的解
% 输入:各步长input_h [hx hy ht theta],差分选项option,是否绘图fig_flag,是否输出误差err_flag,所需误差类型err_types(用""包裹的向量)
% 输出:整个时空上的有限元节点处的解析解,所选err_types的误差err
function [U,err] = FEM_2D_Linear_Dirichlet(input_h,option,fig_flag,err_flag,err_types)
%% 输入网格信息
left = 0;
right = 2;
down = 0;
up = 1;
t_start = 0;
t_end = 1;
boundary_type = [-1 -1 -1 -1];% 边界条件 -1 -> 'D', -2 -> 'Neu', -3 -> 'Robin'

hx = input_h(1);% x 网格步长
hy = input_h(2);% y 网格步长
ht = input_h(3);% t 网格步长
theta = input_h(4);% theta

dim = 2;% 2D
type = 1;% Linear

Time = t_start:ht:t_end;% 时间坐标
[Nf,Nm,Nb,P,T,Pb,Tb,B_m_edges,B_f_edges,f_X,f_Y] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim,type);

%% 组装刚度矩阵与荷载矩阵
M = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @(x,y,t)1, 0, 0, 0, 0, t_start);
A1 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_c, 1, 0, 1, 0, t_start);
A2 = assembler_A(2, Nf, Nf, Nm, P, T, Nb, Nb, [dim type], [dim type], Tb, Tb, @func_c, 0, 1, 0, 1, t_start);
b_all = cell(1,length(Time));
parfor t_ii = 1:length(Time)
    b_all{t_ii} = assembler_b(2, Nf, Nm, P, T, Nb, [dim type], Tb, @func_f, 0, 0, Time(t_ii));
end
A = A1 + A2;
%% 差分求解含时问题
[U,uh_end] = t_fdm(M,A,b_all,theta,Time,func_real_u(Pb(1,:),Pb(2,:),t_start,0,0)',option,P,T,Pb,Tb,B_f_edges,B_m_edges,Nb,Nb,@func_boundary,@func_c,[],[],[], dim, type);

%% 误差分析及可视化
err = [];
if fig_flag == true
    show_eval(option, @func_real_u, U, f_X, f_Y, Time, "Linear 2D FEM")% 调用可视化
end
if err_flag == true
    err = err_eval(@func_real_u, left, right, down, up, uh_end, Nm, Nb, P, T, Pb, Tb, dim, type, input_h(1:3), err_types,t_end);
end
end