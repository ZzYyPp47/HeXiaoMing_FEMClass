function result = basis_ref_function(x, y, dim, type, index, dx, dy)
%% dim = 1 type = 1 1D_Linear FE 
if dim == 1
    if type == 1
        if index == 1
            if dx == 0 
                result = 1 - x;
            elseif dx == 1 
                result = - 1;
            elseif dx >= 2
                result = 0;
            end
        end
        if index == 2
            if dx == 0 
                result = x;
            elseif dx == 1 
                result = 1;
            elseif dx >= 2
                result = 0;
            end
        end
    end
end

%% dim = 1 type = 2 1D_Quadratic FE 
if dim == 1
    if type == 2
        if index == 1
            if dx == 0 
                result = 2.*x.^2-3.*x+1;
            elseif dx == 1 
                result = 4 .* x - 3;
            elseif dx == 2
                result = 4;
            elseif dx >= 3
                result = 0;
            end
        end
        if index == 2
            if dx == 0 
                result = 4 .* x - 4 .* x .^ 2;
            elseif dx == 1 
                result = 4 - 8 .* x;
            elseif dx == 2
                result = -8;
            elseif dx >= 3
                result = 0;
            end
        end
        if index == 3
            if dx == 0 
                result = 2 .* x .^ 2 - x;
            elseif dx == 1 
                result = 4 .* x - 1;
            elseif dx == 2
                result = 4;
            elseif dx >= 3
                result = 0;
            end
        end
    end
end

%% dim = 1 type = 3 1D_Cubic FE 
if dim == 1
    if type == 3
        if index == 1
            if dx == 0 
                result = - (9 .* x .^ 3) ./ 2 + 9 .* x .^ 2 - (11 .* x) ./ 2 + 1;
            elseif dx == 1 
                result = - 27 .* x .^ 2 ./ 2 + 18 .* x - 11 ./ 2;
            elseif dx == 2
                result = -27 .* x + 18;
            elseif dx == 3
                result = -27;
            elseif dx >= 4
                result = 0;
            end
        end
        if index == 2
            if dx == 0 
                result = (27 .* x .^ 3) ./ 2 - (45 .* x .^ 2) ./ 2 + 9 .* x;
            elseif dx == 1 
                result = 81 .* x .^ 2 ./ 2 - 45 .* x + 9;
            elseif dx == 2
                result = 81 .* x - 45;
            elseif dx == 3
                result = 81;
            elseif dx >= 4
                result = 0;
            end
        end
        if index == 3
            if dx == 0 
                result = - (27 .* x .^ 3) ./ 2 + 18 .* x .^ 2 - (9 .* x) ./ 2;
            elseif dx == 1 
                result = -81 .* x .^ 2 ./ 2 + 36 .* x - 9 ./ 2;
            elseif dx == 2
                result = -81 .* x + 36;
            elseif dx == 3
                result = -81;
            elseif dx >= 4
                result = 0;
            end
        end
        if index == 4
            if dx == 0 
                result = 9 .* x .^ 3 ./ 2 - 9 .* x .^ 2 ./ 2 + x;
            elseif dx == 1 
                result = 27 .* x .^ 2 ./ 2 - 9 .* x + 1;
            elseif dx == 2
                result = 27 .* x - 9;
            elseif dx == 3
                result = 27;
            elseif dx >= 4
                result = 0;
            end
        end
    end
end

%% dim = 2 type = 1 2D_Linear FE (triangle)
if dim == 2
    if type == 1
        if index == 1
            if dx == 0 && dy == 0
                result = 1 - x - y;
            elseif dx == 1 && dy == 0
                result = - 1;
            elseif dx == 0 && dy == 1
                result = - 1;
            else
                result = 0;
            end
        end
        if index == 2
            if dx == 0 && dy == 0
                result = x;
            elseif dx == 1 && dy == 0
                result = 1;
            else
                result = 0;
            end
        end
        if index == 3
            if dx == 0 && dy == 0
                result = y;
            elseif dx == 0 && dy == 1
                result = 1;
            else
                result = 0;
            end
        end
    end
end

%% dim = 2 type = 2 2D_Quadratic FE (triangle)
if dim == 2
    if type == 2
        if index == 1
            if dx == 0 && dy == 0
                result = 2 .* x .^ 2 + 2 .* y .^ 2 + 4 .* x .* y - 3 .* y - 3 .* x + 1;
            elseif dx == 1 && dy == 0
                result = 4 .* x + 4 .* y - 3;
            elseif dx == 0 && dy == 1
                result = 4 .* y + 4 .* x - 3;
            elseif dx == 1 && dy == 1
                result = 4;
            else
                result = 0;
            end
        end
        if index == 2
            if dx == 0 && dy == 0
                result = 2 .* x .^ 2 - x;
            elseif dx == 1 && dy == 0
                result = 4 .* x - 1;
            else
                result = 0;
            end
        end
        if index == 3
            if dx == 0 && dy == 0
                result = 2 .* y .^ 2 - y;
            elseif dx == 1 && dy == 0
                result = 0;
            elseif dx == 0 && dy == 1
                result = 4 .* y - 1;
            else
                result = 0;
            end
        end
        if index == 4
            if dx == 0 && dy == 0
                result = -4 .* x .^ 2 - 4 .* x .* y + 4 .* x;
            elseif dx == 1 && dy == 0
                result = -8 .* x - 4 .* y + 4;
            elseif dx == 0 && dy == 1
                result = -4 .* x;
            elseif dx == 1 && dy == 1
                result = -4;
            else
                result = 0;
            end
        end
        if index == 5
            if dx == 0 && dy == 0
                result = 4 .* x .* y;
            elseif dx == 1 && dy == 0
                result = 4 .* y;
            elseif dx == 0 && dy == 1
                result = 4 .* x;
            elseif dx == 1 && dy == 1
                result = 4;
            else
                result = 0;
            end
        end
        if index == 6
            if dx == 0 && dy == 0
                result = -4 .* y .^ 2 - 4 .* x .* y + 4 .* y;
            elseif dx == 1 && dy == 0
                result = -4 .* y;
            elseif dx == 0 && dy == 1
                result = -8 .* y - 4 .* x + 4;
            elseif dx == 1 && dy == 1
                result = -4;
            else
                result = 0;
            end
        end        
    end
end
end