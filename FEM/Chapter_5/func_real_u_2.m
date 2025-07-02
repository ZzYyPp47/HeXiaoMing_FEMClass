function result = func_real_u_2(x,y,dx,dy)
if dx == 0 && dy == 0
    result = x .* (x - 1) .* y .* (y - 1);
elseif dx == 1 && dy == 0
    result = y .* (2 .* x - 1) .* (y - 1);
elseif dx == 0 && dy == 1
    result = x .* (2 .* y - 1) .* (x - 1);
elseif dx == 1 && dy == 1
    result = (2 .* x - 1) .* (2 .* y - 1);
else
    error('Not define dx = %d and dy = %d',dx,dy);
end
end