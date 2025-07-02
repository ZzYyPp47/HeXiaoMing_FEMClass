function result = func_real_u_1(x,y,dx,dy)
if dx == 0 && dy == 0
    result = sin(pi .* x) .* sin(pi .* y);
elseif dx == 1 && dy == 0
    result = pi .* cos(pi .* x) .* sin(pi .* y);
elseif dx == 0 && dy == 1
    result = pi .* cos(pi .* y) .* sin(pi .* x);
elseif dx == 1 && dy == 1
    result = pi .^ 2 .* cos(pi .* x) .* cos(pi .* y);
else
    error('Not define dx = %d and dy = %d',dx,dy);
end
end