function result = func_f2(x,y,t)
result = (4.*func_nu(x,y).*x.*y-func_nu(x,y).*pi.^3.*sin(pi.*x)+2.*pi.*(2-pi.*sin(pi.*x)).*sin(2.*pi.*y)).*cos(2.*pi.*t)-2.*pi.*(-2.*x.*y.^3./3+2-pi.*sin(pi.*x)).*sin(2.*pi.*t) + ((x.^2.*y.^2+exp(-y)).*(-2.*y.^3./3-pi.^2.*cos(pi.*x))+(-2.*x.*y.^3./3+2-pi.*sin(pi.*x)).*(-2.*x.*y.^2)).*cos(2.*pi.*t).*cos(2.*pi.*t);
end