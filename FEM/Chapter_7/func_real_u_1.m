function result = func_real_u_1(x,y,dx,dy)
if dx == 0 && dy == 0
    result = x .^ 2 .* y .^ 2 + exp(-y);
elseif dx == 1 && dy == 0
    result = 2 .* x .* y .^ 2;
elseif dx == 0 && dy == 1
    result = 2 .* y .* x .^ 2 - exp(-y);
elseif dx == 1 && dy == 1
    result = 4 .* x .* y;
else
    error('Not define dx = %d and dy = %d',dx,dy);
end
end