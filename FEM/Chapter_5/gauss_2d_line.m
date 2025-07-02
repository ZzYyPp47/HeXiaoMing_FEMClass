% 高斯-勒让德求积(optional) 2D 第一型曲线积分
% f is a function handle which can feed in (x,y)
function result = gauss_2d_line(f,a,b,P_1,P_2,type)
gauss_points = [0,-0.3242534234,-0.6133714327,-0.8360311073,-0.9681602395,0.3242534234,0.6133714327, 0.8360311073, 0.9681602395];
gauss_weights = [0.330239355,0.312347077,0.2606106964,0.1806481607,0.08127438836,0.312347077,0.2606106964,0.1806481607,0.08127438836];
if type == 'v'% 积分线垂直于x轴,此时a b为[ymin ymax],line_func为积分曲线@(x)c0
    gauss_mapped_pointed = (b - a) .* gauss_points ./ 2 + (a + b) ./ 2;
    result = (b - a)/2 * gauss_weights * f(line_func(gauss_mapped_pointed, P_1, P_2),gauss_mapped_pointed)';
else% 其他情况,此时a b为[xmin xmax],line_func为积分曲线@(x)ax+b
    gauss_mapped_pointed = (b - a) .* gauss_points ./ 2 + (a + b) ./ 2;
    k = (P_1(2) - P_2(2)) / (P_1(1) - P_2(1));
    result = sqrt(1 + k^2) * (b - a)/2 * gauss_weights * f(gauss_mapped_pointed,line_func(gauss_mapped_pointed, P_1, P_2))';
end