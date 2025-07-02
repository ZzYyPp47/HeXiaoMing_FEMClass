function result = func_boundary(x, y)
% point = [x y]
% 支持向量化输入

% 初始化结果数组
result = zeros(size(x));

% 逻辑判断并赋值
result(x == -1) = exp(-1 + y(x == -1));
result(x == 1) = exp(1 + y(x == 1));
% result(y == -1) = -2 * x(y == -1) .* (1 - x(y == -1) / 2) .* exp(x(y == -1) - 1);
result(y == 1) = exp(x(y == 1) + 1);

% 检查是否有不在边界上的点
if any(~(x == -1 | x == 1 | y == -1 | y == 1))
    error('不在边界上.');
end
end