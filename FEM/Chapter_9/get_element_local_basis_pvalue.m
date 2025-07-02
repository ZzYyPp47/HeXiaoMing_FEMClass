% use to evaluate arbitrary point in one element
% WILL NOT CHECK WHETHER THE POINT IN ELEMENT OR NOT.
function result = get_element_local_basis_pvalue(uh, x, y, element_index, Nb, P, T, Tb, dim, type, dx, dy)
uh_local = uh(Tb(:,element_index));
temp = zeros([size(x),Nb]);% 存储评估点在各局部基函数上的值
for ii = 1:1:Nb% 遍历所有局部基函数
    temp(:,:,ii) = uh_local(ii) .* basis_local_function(x, y, dim, type, ii, dx, dy, P(:,T(:,element_index)));
end
result = sum(temp,3);
end