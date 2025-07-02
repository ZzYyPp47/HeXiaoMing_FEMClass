function result = func_real_p(x,y,t,dx,dy)
if dx == 0 && dy == 0
    result = cos(2.*pi.*t).*cos(2.*pi.*y).*(pi.*sin(pi.*x) - 2);
elseif dx == 1 && dy == 0
    result = pi.^2.*cos(2.*pi.*t).*cos(pi.*x).*cos(2.*pi.*y);
elseif dx == 0 && dy == 1
    result = -2.*pi.*cos(2.*pi.*t).*sin(2.*pi.*y).*(pi.*sin(pi.*x) - 2);
elseif dx == 1 && dy == 1
    result = -2.*pi.^3.*cos(2.*pi.*t).*cos(pi.*x).*sin(2.*pi.*y);
else
    error('Not define dx = %d and dy = %d',dx,dy);
end
end