function result = func_boundary2(x, y)
% point = [x y]
% 支持向量化输入

% 初始化结果数组
result = zeros(size(x));

% 逻辑判断并赋值
result(x == 0) = 2;
result(x == 1) = -2 .* y(x == 1) .^ 3 ./ 3 + 2;
result(y == -0.25) = x(y == -0.25) ./ 96 + 2 - pi .* sin(pi .* x(y == -0.25));
result(y == 0) = 2 - pi .* sin(pi .* x(y == 0));

% 检查是否有不在边界上的点
if any(~(x == 0 | x == 1 | y == 0 | y == -0.25))
    error('不在边界上.');
end
end