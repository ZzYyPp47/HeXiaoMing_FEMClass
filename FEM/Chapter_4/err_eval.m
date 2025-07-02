% max,L^infty,L^2,H^1(semi),H^1误差
% h = [hx hy]
function err = err_eval(real_u, left, right, down, up, uh, Nm, Nb, P, T, Pb, Tb, dim, type, h, err_types, t0)
err = h;
F = @(x,y,element_index,dx,dy)(func_real_u(x,y,t0,dx,dy) - get_element_local_basis_pvalue(uh, x, y, element_index, Nb, P, T, Tb, dim, type, dx, dy)) .^2;
G = @(x,y,element_index)max(abs(func_real_u(x,y,t0,0,0) - get_element_local_basis_pvalue(uh, x, y, element_index, Nb, P, T, Tb, dim, type, 0, 0)));
for ii = 1:1:length(err_types)
    switch err_types(ii)
        case "max"% 节点max误差
            temp = max(abs(real_u(Pb(1,:),Pb(2,:),t0,0,0)' - uh));
        case "L^inf"% L^infty误差
            temp_all = zeros(1,Nm);
            parfor element_index = 1:1:Nm
                eval_points = get_element_gauss_point(2, P(:,T(:,element_index)));
                temp_all(element_index) = G(eval_points(1,:),eval_points(2,:),element_index);
            end
            temp = max(temp_all);
        case "L^2"% L^2误差
            temp = 0;
            parfor element_index = 1:1:Nm
                f_1 = @(x,y)F(x,y,element_index,0,0);
                sum = gauss_2d_triangle(f_1,P(:,T(:,element_index)));
                temp = temp + sum;
            end
            temp = sqrt(temp);
        case "H^1(semi)"% H^1(semi)误差
            temp = 0;
            parfor element_index = 1:1:Nm
                f_1 = @(x,y)F(x,y,element_index,1,0) + F(x,y,element_index,0,1);
                sum = gauss_2d_triangle(f_1,P(:,T(:,element_index)));
                temp = temp + sum;
            end
            temp = sqrt(temp);
        case "H^1"% H^1误差
            temp = 0;
            parfor element_index = 1:1:Nm
                f_1 = @(x,y)F(x,y,element_index,1,0) + F(x,y,element_index,0,1) + F(x,y,element_index,0,0);
                sum = gauss_2d_triangle(f_1,P(:,T(:,element_index)));
                temp = temp + sum;
            end
            temp = sqrt(temp);
        case "all"
            temp = err_eval(real_u, left, right, down, up, uh, Nm, Nb, P, T, Pb, Tb, dim, type, h, ["max","L^inf","L^2","H^1","H^1(semi)"],t0);
            temp = temp(4:end);% 去掉重复的 h 
        otherwise
            error(['Unsupported error type: ' err_types(ii)]);
    end
    err = [err temp];
end
end