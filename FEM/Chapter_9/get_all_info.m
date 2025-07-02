function [Nf,Nm,Nb,P,T,Pb,Tb,B_m_edges,B_f_edges,f_X,f_Y] = get_all_info(left,right,down,up,hx,hy,boundary_type,dim,type)
%% 获取网格信息
Nx = (right - left) / hx;% x网格单元数
Ny = (up - down) / hy;% y网格单元数

%% 2D_Linear
if dim == 2 && type == 1
% 获取有限元信息
Nfx = Nx + 1;% x有限元节点数
Nfy = Ny + 1;% y有限元节点数
Nf = Nfx * Nfy;% 有限元节点数
Nm = 2 * Nx * Ny;% 网格单元数
Nb = 3;% 局部基函数个数

% 初始化其余信息
[m_idx_X,m_idx_Y] = meshgrid(1:1:Nx,1:1:Ny);% 网格单元自然坐标
[m_X,m_Y] = meshgrid(left:hx:right,down:hy:up);% 网格坐标
[f_idx_X,f_idx_Y] = meshgrid(1:1:Nfx,1:1:Nfy);% 有限元节点自然坐标
[f_X,f_Y] = deal(m_X,m_Y);% 有限元坐标

% 构建P,T,Pb,Tb信息矩阵
P = [reshape(m_X,1,[]);reshape(m_Y,1,[])];% 网格点坐标矩阵
T = m_index_to_tri_T(reshape(m_idx_X,1,[]),reshape(m_idx_Y,1,[]),Ny);% 网格点索引矩阵
Pb = P;% 有限元点坐标矩阵
Tb = T;% 有限元点索引矩阵
B_m_edges = edges_to_tri_mbidx(m_idx_X,m_idx_Y,T,Nx,Ny,boundary_type);% 网格边界信息矩阵(第一行是边界类型,第二行是网格边索引,其余行是网格边端点)
B_f_edges = edges_to_tri_fbidx(f_idx_X,f_idx_Y,Nfx,Nfy,boundary_type);% 有限元边界信息矩阵(第一行是边界类型,第二行是有限元边节点)
B_f_edges = check_legal(B_m_edges,B_f_edges,Nx,Ny,Nfx,Nfy);% 确保正确施加边界条件
end

%% 2D_Quadratic
if dim == 2 && type == 2
% 获取有限元信息
Nfx = 2 * Nx + 1;% x有限元节点数
Nfy = 2 * Ny + 1;% y有限元节点数
Nf = Nfx * Nfy;% 有限元节点数
Nm = 2 * Nx * Ny;% 网格单元数
Nb = 6;% 局部基函数个数

% 初始化其余信息
[m_idx_X,m_idx_Y] = meshgrid(1:1:Nx,1:1:Ny);% 网格单元自然坐标
[m_X,m_Y] = meshgrid(left:hx:right,down:hy:up);% 网格坐标
[f_idx_X,f_idx_Y] = meshgrid(1:1:Nfx,1:1:Nfy);% 有限元节点自然坐标
[f_X,f_Y] = meshgrid(left:hx/2:right,down:hy/2:up);% 有限元坐标

% 构建P,T,Pb,Tb信息矩阵
P = [reshape(m_X,1,[]);reshape(m_Y,1,[])];% 网格点坐标矩阵
T = m_index_to_tri_T(reshape(m_idx_X,1,[]),reshape(m_idx_Y,1,[]),Ny);% 网格点索引矩阵
Pb = [reshape(f_X,1,[]);reshape(f_Y,1,[])];% 有限元点坐标矩阵
Tb = f_index_to_qua_Tb(reshape(m_idx_X,1,[]),reshape(m_idx_Y,1,[]),Nfy,Ny);% 有限元点索引矩阵
B_m_edges = edges_to_tri_mbidx(m_idx_X,m_idx_Y,T,Nx,Ny,boundary_type);% 网格边界信息矩阵(第一行是边界类型,第二行是网格边索引,其余行是网格边端点)
B_f_edges = edges_to_tri_fbidx(f_idx_X,f_idx_Y,Nfx,Nfy,boundary_type);% 有限元边界信息矩阵(第一行是边界类型,第二行是有限元边节点)
B_f_edges = check_legal(B_m_edges,B_f_edges,Nx,Ny,Nfx,Nfy);% 确保正确施加边界条件
end

