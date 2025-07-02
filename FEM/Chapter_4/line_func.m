% Get a line function from given two points
% P_1 = [x1 y1]; P_2 = [x2 y2];
function result = line_func(x, P_1, P_2)
if P_1(1) == P_2(1)
    result = P_1(1);
elseif P_1(2) == P_2(2)
    result = P_2(2);
else
    result = (P_1(2) - P_2(2)) / (P_1(1) - P_2(1)) * (x - P_1(1)) + P_1(2);
end