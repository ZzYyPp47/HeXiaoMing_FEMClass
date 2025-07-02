function result = func_f1(x,y)
result = -(func_lamda(x,y) + 3 .* func_nu(x,y)) .* (-pi^2 .* sin(pi .* x) .* sin(pi .* y))...
         -(func_lamda(x,y) + func_nu(x,y)) .* ((2 .* x - 1) .* (2 .* y - 1));
end