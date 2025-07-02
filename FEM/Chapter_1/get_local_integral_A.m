% compute each integral
% basis_trial/test = [dim type]
function result = get_local_integral_A(int_type, alpha, beta, basis_trial, basis_test, func_c, trial_dx, trial_dy, test_dx, test_dy, vertices)

% do 1D integral, vertices = [x1,x2]
if int_type == 1 
    int_func = @(x) func_c(x, [])...
                 .* basis_local_function(x, [], basis_trial(1), basis_trial(2), alpha, trial_dx, trial_dy, vertices)...
                 .* basis_local_function(x, [], basis_test(1), basis_test(2), beta, test_dx, test_dy, vertices);
    result = gauss_1d(int_func, vertices);
end

% do 2D integral (triangle), vertices = [x1 x2 x3;y1 y2 y3]
if int_type == 2 
    int_func = @(x,y) func_c(x, y)...
                 .* basis_local_function(x, y, basis_trial(1), basis_trial(2), alpha, trial_dx, trial_dy, vertices)...
                 .* basis_local_function(x, y, basis_test(1), basis_test(2), beta, test_dx, test_dy, vertices);
    result = gauss_2d_triangle(int_func, vertices);
end
end

