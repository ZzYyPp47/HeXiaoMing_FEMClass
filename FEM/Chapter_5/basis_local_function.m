function result = basis_local_function(x, y, dim, type, index, dx, dy, vertices)
%% dim = 1  1D_Linear FE && 1D_Quadratic FE
% vertices = [x1,x2]
if dim == 1
    h = max(vertices) - min(vertices);
    x_hat = (x - min(vertices)) ./ h;
    if dx == 0 
        result = basis_ref_function(x_hat, y, dim, type, index, dx, dy);
    elseif dx >= 1 
        result = basis_ref_function(x_hat, y, dim, type, index, dx, dy) ./ h;% d2x_hat/dx_hat^2 = 0
    end
end

%% dim = 2 2D_Linear FE && 2D_Quadratic FE (triangle)
% vertices = [x1 x2 x3; y1 y2 y3]
if dim == 2
    det_J = det([vertices(:,2) - vertices(:,1),vertices(:,3) - vertices(:,1)]);
    x_hat = ((vertices(2,3) - vertices(2,1)) .* (x - vertices(1,1)) - (vertices(1,3) - vertices(1,1)) .* (y - vertices(2,1))) ./ det_J;
    y_hat = ((vertices(1,2) - vertices(1,1)) .* (y - vertices(2,1)) - (vertices(2,2) - vertices(2,1)) .* (x - vertices(1,1))) ./ det_J;
    x_13 = (vertices(1,1) - vertices(1,3)) ./ det_J;
    x_21 = (vertices(1,2) - vertices(1,1)) ./ det_J;
    y_31 = (vertices(2,3) - vertices(2,1)) ./ det_J; 
    y_12 = (vertices(2,1) - vertices(2,2)) ./ det_J;
    if dx == 0 && dy == 0
        result = basis_ref_function(x_hat, y_hat, dim, type, index, dx, dy);
    elseif dx == 1 && dy == 0
        result = basis_ref_function(x_hat, y_hat, dim, type, index, 1, 0) .* y_31...
               + basis_ref_function(x_hat, y_hat, dim, type, index, 0, 1) .* y_12;
    elseif dx == 0 && dy == 1
        result = basis_ref_function(x_hat, y_hat, dim, type, index, 1, 0) .* x_13...
               + basis_ref_function(x_hat, y_hat, dim, type, index, 0, 1) .* x_21;
    elseif dx == 1 && dy == 1
        result = basis_ref_function(x_hat, y_hat, dim, type, index, 2, 0) .* x_13 .* y_31...
               + basis_ref_function(x_hat, y_hat, dim, type, index, 1, 1) .* x_13 .* y_12...
               + basis_ref_function(x_hat, y_hat, dim, type, index, 1, 1) .* x_21 .* y_31...
               + basis_ref_function(x_hat, y_hat, dim, type, index, 0, 2) .* x_21 .* y_12;
    elseif dx == 2 && dy == 0
        result = basis_ref_function(x_hat, y_hat, dim, type, index, 2, 0) .* y_31 .* y_31...
               + 2 .* basis_ref_function(x_hat, y_hat, dim, type, index, 1, 1) .* y_31 .* y_12...
               + basis_ref_function(x_hat, y_hat, dim, type, index, 0, 2) .* y_12 .* y_12;
    elseif dx == 0 && dy == 2
        result = basis_ref_function(x_hat, y_hat, dim, type, index, 2, 0) .* x_13 .* x_13...
               + 2 .* basis_ref_function(x_hat, y_hat, dim, type, index, 1, 1) .* x_13 .* x_21...
               + basis_ref_function(x_hat, y_hat, dim, type, index, 0, 2) .* x_21 .* x_21;
    else
        result = 0;
    end
end
end