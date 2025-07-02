% vertices = [x1 x2] in 1D or vertices = [x1 x2 x3;y1 y2 y3] in 2D
function mapped_gauss_points = get_element_gauss_point(dim, vertices)
gauss_reference_points_1d = [0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425];
gauss_reference_points_2d = [1/2,(1+sqrt(3/5))/2,(1+sqrt(3/5))/2,(1-sqrt(3/5))/2,(1-sqrt(3/5))/2, 1/2, 1/2, (1+sqrt(3/5))/2,(1-sqrt(3/5))/2;...
    1/4,(1-sqrt(3/5))*(1+sqrt(3/5))/4,(1-sqrt(3/5))*(1-sqrt(3/5))/4,(1+sqrt(3/5))*(1+sqrt(3/5))/4,(1+sqrt(3/5))*(1-sqrt(3/5))/4,(1+sqrt(3/5))/4,(1-sqrt(3/5))/4,(1-sqrt(3/5))/4,(1+sqrt(3/5))/4];
if dim == 1
    [xmax,xmin] = deal(max(vertices),min(vertices));
    mapped_gauss_points = (xmax - xmin) ./ 2 .* gauss_reference_points_1d + (xmax + xmin) ./ 2;
elseif dim == 2
    gauss_points_x = vertices(1,1) .* shape_func(gauss_reference_points_2d(1,:),gauss_reference_points_2d(2,:),1)...
                   + vertices(1,2) .* shape_func(gauss_reference_points_2d(1,:),gauss_reference_points_2d(2,:),2)...
                   + vertices(1,3) .* shape_func(gauss_reference_points_2d(1,:),gauss_reference_points_2d(2,:),3);
    gauss_points_y = vertices(2,1) .* shape_func(gauss_reference_points_2d(1,:),gauss_reference_points_2d(2,:),1)...
                   + vertices(2,2) .* shape_func(gauss_reference_points_2d(1,:),gauss_reference_points_2d(2,:),2)...
                   + vertices(2,3) .* shape_func(gauss_reference_points_2d(1,:),gauss_reference_points_2d(2,:),3);
    mapped_gauss_points = [gauss_points_x;gauss_points_y];
end
end

function result = shape_func(x,y,type)
if type == 1
    result = 1 - x - y;
elseif type == 2
    result = x;
elseif type == 3
    result = y;
end
end