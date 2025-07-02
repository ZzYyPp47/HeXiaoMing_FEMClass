% point = [x1 x2 ...] in 1D or point = [x1 x2 x3...;y1 y2 y3...] in 2D
function result = get_uh_point(uh, point, Nm, Nb, P, T, Tb, dim, type, dx, dy)
result = zeros(1,size(point,2));% 存储评估点处的值
for ii = 1:1:size(point,2)
    for jj = 1:1:Nm% 跑遍每个单元
        if dim == 1
            if point(ii) >= P(jj) && point(ii) <= P(jj + 1) % 在 E_jj 上
                result(ii) = get_element_local_basis_pvalue(uh, point(ii), [], jj, Nb, P, T, Tb, dim, type, dx, dy);
                break;% 继续下个评估点
            end
        else % dim == 2
            if inpolygon(point(1,ii), point(2,ii), P(1,T(:,jj)), P(2,T(:,jj))) % 判断是否在 E_jj 上
                result(ii) = get_element_local_basis_pvalue(uh, point(1,ii), point(2,ii), jj, Nb, P, T, Tb, dim, type, dx, dy);
                break; % 继续下个评估点
            end
        end
    end
end
end