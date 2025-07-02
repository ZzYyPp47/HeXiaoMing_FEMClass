function result = func_real_u_2(x,y,t,dx,dy)
if dx == 0 && dy == 0
    result = (2 - (2 .* x .* y .^ 3) ./ 3 - pi .* sin(pi .* x)).*cos(2.*pi.*t);
elseif dx == 1 && dy == 0
    result = (- pi .^2 .* cos(pi .* x) - (2 .* y .^ 3) ./ 3).*cos(2.*pi.*t);
elseif dx == 0 && dy == 1
    result = (-2 .* x .* y .^ 2).*cos(2.*pi.*t);
elseif dx == 1 && dy == 1
    result = (-2 .* y .^ 2).*cos(2.*pi.*t);
else
    error('Not define dx = %d and dy = %d',dx,dy);
end
end