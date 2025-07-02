% Convert finite nodes natural index to boundary 1D finite nodes index with boundary type
% Local index: down -> right -> up -> left
function B = edges_to_tri_fbidx(f_idx_X,f_idx_Y,Nfx,Nfy,boundary_type)
to_idx = @(x,y)(x - 1) * Nfy + y;% Convert finite nodes natural index to 1D index
down_f_idx = to_idx(reshape(f_idx_X(1,:),1,[]),reshape(f_idx_Y(1,:),1,[]));
right_f_idx = to_idx(reshape(f_idx_X(2:end,end),1,[]),reshape(f_idx_Y(2:end,end),1,[]));
up_f_idx = flip(to_idx(reshape(f_idx_X(end,1:end - 1),1,[]),reshape(f_idx_Y(end,1:end - 1),1,[])));
left_f_idx = flip(to_idx(reshape(f_idx_X(2:end - 1,1),1,[]),reshape(f_idx_Y(2:end - 1,1),1,[])));
temp = [down_f_idx,right_f_idx,up_f_idx,left_f_idx];
B = [repelem(boundary_type,[Nfx,Nfy - 1,Nfx - 1,Nfy - 2]);temp];
end