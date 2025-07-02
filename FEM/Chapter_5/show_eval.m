% 可视化函数
% compute finite nodes error.
% so we feed in point_x point_y as finite coordinate.
function show_eval(real_u, uh, Pb_trial, point_x, point_y, txt)
real_y = reshape(real_u(Pb_trial(1,:),Pb_trial(2,:),0,0),size(point_x));
uh_y = reshape(uh,size(point_x));
level = 100;

% 比较
figure;
tiledlayout(2,1);
nexttile;
contourf(point_x,point_y,real_y,level,'EdgeAlpha',0)
title('Real');
xlabel('x');
ylabel('y');
cb = colorbar;
cb.Layout.Tile = 'east';

nexttile;
contourf(point_x,point_y,uh_y,level,'EdgeAlpha',0)
title('FEM');
xlabel('x');
ylabel('y');

% 误差分析
abs_error = abs(real_y - uh_y);
figure;
contourf(point_x,point_y,abs_error,level,'EdgeAlpha',0)
title("Absolute Error (" + txt + ")");
xlabel('x');
ylabel('y');
colorbar;
end