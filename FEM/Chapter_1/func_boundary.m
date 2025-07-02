function result = func_boundary(x, y)
% point = [x y]
% 支持向量化输入

% 初始化结果数组
result = zeros(size(x));

% 逻辑判断并赋值
result(x == 0) = 0;
result(x == 1) = cos(1);

% 检查是否有不在边界上的点
if any(~(x == 0 | x == 1))
    error('不在边界上.');
end
end