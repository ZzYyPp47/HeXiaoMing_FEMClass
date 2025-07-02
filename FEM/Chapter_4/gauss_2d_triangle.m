% func need to be a function handle which can feed in (x,y).
% vertices is a form like [x1 x2 x3...;y1 y2 y3...]
function result = gauss_2d_triangle(func,vertices)
gauss_reference_weights = [64/81*(1-0)/8,100/324*(1-sqrt(3/5))/8,100/324*(1-sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,100/324*(1+sqrt(3/5))/8,40/81*(1-0)/8,40/81*(1-0)/8,40/81*(1-sqrt(3/5))/8,40/81*(1+sqrt(3/5))/8];
gauss_reference_points = [1/2,(1+sqrt(3/5))/2,(1+sqrt(3/5))/2,(1-sqrt(3/5))/2,(1-sqrt(3/5))/2, 1/2, 1/2, (1+sqrt(3/5))/2,(1-sqrt(3/5))/2;...
                          1/4,(1-sqrt(3/5))*(1+sqrt(3/5))/4,(1-sqrt(3/5))*(1-sqrt(3/5))/4,(1+sqrt(3/5))*(1+sqrt(3/5))/4,(1+sqrt(3/5))*(1-sqrt(3/5))/4,(1+sqrt(3/5))/4,(1-sqrt(3/5))/4,(1-sqrt(3/5))/4,(1+sqrt(3/5))/4];
det_J = det([vertices(:,2) - vertices(:,1),vertices(:,3) - vertices(:,1)]);
gauss_points_x = vertices(1,1) .* shape_func(gauss_reference_points(1,:),gauss_reference_points(2,:),1)...
               + vertices(1,2) .* shape_func(gauss_reference_points(1,:),gauss_reference_points(2,:),2)...
               + vertices(1,3) .* shape_func(gauss_reference_points(1,:),gauss_reference_points(2,:),3);
gauss_points_y = vertices(2,1) .* shape_func(gauss_reference_points(1,:),gauss_reference_points(2,:),1)...
               + vertices(2,2) .* shape_func(gauss_reference_points(1,:),gauss_reference_points(2,:),2)...
               + vertices(2,3) .* shape_func(gauss_reference_points(1,:),gauss_reference_points(2,:),3);   
result = sum(gauss_reference_weights .* det_J .* func(gauss_points_x,gauss_points_y),"all");
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
