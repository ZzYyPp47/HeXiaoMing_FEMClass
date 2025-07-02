% n 为求积阶数
function [point,weight] = get_gauss_point_weight(n)
syms x;
N = 10;% 有效数字位数
Poly = diff((x.^2 - 1).^(n + 1),x,n + 1) ./ (factorial(n + 1) .* 2 .^ (n + 1));
Poly_d = matlabFunction(diff(Poly,x));
point = vpa(solve(Poly),N);
weight = vpa(2 ./ ((1 - point.^2).*(Poly_d(point).^2)),N);
end