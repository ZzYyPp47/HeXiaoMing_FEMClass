% Convert T to boundary edges mesh nodes with boundary type
% Local index: down -> right -> up -> left
function B = edges_to_tri_mbidx(e_idx_X,e_idx_Y,T,Nx,Ny,boundary_type)
to_idx_l = @(x,y)(x - 1) * 2 * Ny + 2 * y - 1;% Convert natural elements index to 1D index (left)
to_idx_r = @(x,y)(x - 1) * 2 * Ny + 2 * y;% Convert natural elements index to 1D index (right)
down_m_idx = to_idx_l(reshape(e_idx_X(1,:),1,[]),reshape(e_idx_Y(1,:),1,[]));
right_m_idx = to_idx_r(reshape(e_idx_X(:,end),1,[]),reshape(e_idx_Y(:,end),1,[]));
up_m_idx = flip(to_idx_r(reshape(e_idx_X(end,:),1,[]),reshape(e_idx_Y(end,:),1,[])));
left_m_idx = flip(to_idx_l(reshape(e_idx_X(:,1),1,[]),reshape(e_idx_Y(:,1),1,[])));
temp = [down_m_idx,right_m_idx,up_m_idx,left_m_idx;...
        T([1 2],down_m_idx),T([2 3],right_m_idx),T([3 1],up_m_idx),T([3 1],left_m_idx)];% find corresponding mesh nodes 1d index
B = [repelem(boundary_type,[Nx,Ny,Nx,Ny]);temp];
end