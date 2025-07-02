% Convert natural elements index to matrix T
function T = m_index_to_tri_T(flatten_eX,flatten_eY,Ny)
f = @(x,y) (Ny + 1) * (x - 1) + y;% Convert natural mesh nodes index to 1D index
to_T = @(x,y) [(x - 1) * 2 * Ny + 2 * y - 1,(x - 1) * 2 * Ny + 2 * y;...
               f(x,y),f(x,y + 1);...
               f(x + 1,y),f(x + 1,y);...
               f(x,y + 1),f(x + 1,y + 1)];% natural elements index -> 1D elements index & natural mesh nodes index -> 1D mesh nodes index
T = to_T(flatten_eX,flatten_eY);
T(:,T(1,:)) = T(:,sort(T(1,:)));% 按照自然顺序
T = T(2:end,:);% 抹去第一列
end