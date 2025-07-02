function result = func_f2(x,y)
result = 4 .* func_nu(x,y) .* x .* y - func_nu(x,y) .* pi .^ 3 .* sin(pi .* x)...
           + 2 .* pi .* (2 - pi .* sin(pi .* x)) .* sin(2 .* pi .* y);
end