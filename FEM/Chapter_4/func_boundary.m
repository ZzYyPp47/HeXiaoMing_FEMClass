function result = func_boundary(x, y, t0)
% point = [x y t]
% 支持向量化输入

% 初始化结果数组
result = zeros(size(x));

% 逻辑判断并赋值
result(x == 0) = exp(y(x == 0) + t0);
result(x == 2) = exp(2 + y(x == 2) + t0);
result(y == 0) = exp(x(y == 0) + t0);
result(y == 1) = exp(x(y == 1) + 1 + t0);

% 检查是否有不在边界上的点
if any(~(x == 0 | x == 2 | y == 0 | y == 1))
    error('不在边界上.');
end
end