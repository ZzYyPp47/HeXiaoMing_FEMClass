function result = func_real_u(x,y,dx,dy)
if dx == 0 
    result = x .* cos(x);
elseif dx == 1 
    result = cos(x) - x .* sin(x);
else
    error('Not define dx = %d and dy = %d',dx,dy);
end
end