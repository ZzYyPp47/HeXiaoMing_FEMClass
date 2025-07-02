function result = func_f(x,y)
result = -y.*(1-y).*(1-x-x.^2./2).*exp(x+y)-x.*(1-x./2).*(-3.*y-y.^2).*exp(x+y);
end