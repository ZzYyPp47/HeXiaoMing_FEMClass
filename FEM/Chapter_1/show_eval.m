% 可视化函数
% compute finite nodes error.
% so we feed in point_x point_y as finite coordinate.
function show_eval(real_u, uh, Pb, txt)
uh_y = uh;
% 比较
figure;
plot(Pb,real_u(Pb,[],0,0),'LineWidth',2);
hold on;
plot(Pb,uh_y,'--','LineWidth',2);
grid on;
grid minor;
title("Comparsion (" + txt + ")");
xlabel('x');
ylabel('y');
legend('Numerical Solution','Analytic Solution');
% 误差分析
abs_error = abs(real_u(Pb,[],0,0)' - uh_y);
figure;
plot(Pb,abs_error,'LineWidth',2);
grid on;
grid minor;
title("Absolute Error (" + txt + ")");
xlabel('x');
ylabel('Absolute Error');
end