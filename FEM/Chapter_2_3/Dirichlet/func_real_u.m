function result = func_real_u(x,y,dx,dy)
if dx == 0 && dy == 0
    result = x.*y.*(1-x./2).*(1-y).*exp(x+y);
elseif dx == 1 && dy == 0
    result = (y.*exp(x + y).*(x.^2 - 2).*(y - 1))./2;
elseif dx == 0 && dy == 1
    result = (x.*exp(x + y).*(x - 2).*(y.^2 + y - 1))./2;
elseif dx == 1 && dy == 1
    result = (exp(x + y).*(x.^2 - 2).*(y.^2 + y - 1))./2;
else
    error('Not define dx = %d and dy = %d',dx,dy);
end
end