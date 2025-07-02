function result = func_f1(x,y)
result = -2 .* func_nu(x,y) .* x .^ 2 - 2 .* func_nu(x,y) .* y .^ 2 ...
           - func_nu(x,y) .* exp(-y) + pi .^ 2 .* cos(pi .* x) .* cos(2 .* pi .* y);
end