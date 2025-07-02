function result = func_f1(x,y,t)
result = (-2.*func_nu(x,y).*x.^2-2.*func_nu(x,y).*y.^2-func_nu(x,y).*exp(-y)+pi.^2.*cos(pi.*x).*cos(2.*pi.*y)).*cos(2.*pi.*t)-2.*pi.*(x.^2.*y.^2+exp(-y)).*sin(2.*pi.*t)+((x.^2.*y.^2+exp(-y)).*2.*x.*y.^2+(-2.*x.*y.^3./3+2-pi.*sin(pi.*x)).*(2.*x.^2.*y-exp(-y))).*cos(2.*pi.*t).*cos(2.*pi.*t);
end