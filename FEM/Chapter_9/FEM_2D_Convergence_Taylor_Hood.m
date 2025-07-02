%% Covergence Test for 2D example FEM
% 计算 (u1 u2)t + (u1 u2 · ∇)(u1 u2) - ∇· T(u1 u2 p) = (f1 f2) on \Omega x T
% ∇· (u1 u2) = (0 0) on \Omega x T
% u1 = x^2*y^2 + exp(-y), t = 0;
% u2 = -2*x*y^3/3 + 2 - pi * sin(pi * x),-(2 - pi * sin(pi * x)) * cos(2 * pi * y) t = 0;
% p = -(2 - pi * sin(pi * x)) * cos(2 * pi * y) t = 0;
% f1 = -2 .* pi .* (x.^2.*y.^2+exp(-y)).*sin(2.*pi.*t)+(-2 * mu * x^2 - 2 * mu * y ^2 - mu*exp(-y) + pi^2 * cos(pi * x) * cos(2 * pi * y)).*cos(2.*pi.*t);
% f2 = -2.*pi.*(-2.*x.*y.^3./3+2-pi.*sin(pi.*x)).*sin(2.*pi.*t)+(4 * mu * x * y - mu * pi^3 * sin(pi * x) + 2* pi * (2 - pi * sin(pi * x)) * sin(2 * pi * y)).*cos(2.*pi.*t);
% mu = 1;
% real solution: [u1 u2 p] = [(x^2*y^2 + exp(-y)).*cos(2.*pi.*t),(-2*x*y^3/3 + 2 - pi * sin(pi * x)).*cos(2.*pi.*t),-(2 - pi * sin(pi * x)) * cos(2 * pi * y)*cos(2*pi*t)]
% \Omega = [0,1] x [-0.25,0] x [0,1]
%% 初始化
clc;
clear;
N = 3;
result_u = zeros(N,7);
result_p = zeros(N,7);
f_1 = @(x)x;
f_2 = @(x)x.^2;
f_3 = @(x)x.^3;
poolobj = gcp('nocreate');
if isempty(poolobj)
    threads = 32; % 设定机器线程数
    fprintf('Initialize parallel pool first, with threads %d.\n',threads)
    parpool(threads);
    fprintf('Now starting...\n');
else
    fprintf('Already initialized parallel pool, with threads %d.\n',poolobj.NumWorkers);
    fprintf('Now starting...\n');
end
%% 计算
t_s = tic;
for ii = 1:N
    [~,result_u(ii,:),result_p(ii,:)] = FEM_2D_Taylor_Hood_Dirichlet([1/2^(ii + 2),1/2^(ii + 2),8 * (1/2^(ii + 2)) ^ 3],"end",false,true,["max","L^inf","L^2","H^1(semi)"]);
    time = toc(t_s);
    fprintf('Completed task num: %d, time passed: %.8f s, Total: %d\n',ii,time,N);
end
% u
result_order_u = log(result_u(1:end - 1,4:end)./result_u(2:end,4:end)) ./ log(result_u(1:end - 1,1)./result_u(2:end,1));
result_order_u = table(result_order_u(:,1),result_order_u(:,2),result_order_u(:,3),result_order_u(:,4),'VariableNames',{'Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
draw(f_2,f_3,result_u,['(Quadratic 2D FEM)',"h^2","h^3"]);
result_u = table(sym(result_u(:,1)),sym(result_u(:,2)),sym(result_u(:,3)),result_u(:,4),result_u(:,5),result_u(:,6),result_u(:,7),...
               'VariableNames',{'hx','hy','ht','Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
fprintf("u: \n");
disp(result_u);
fprintf("Convergence Order: \n");
disp(result_order_u);
% p
result_order_p = log(result_p(1:end - 1,4:end)./result_p(2:end,4:end)) ./ log(result_p(1:end - 1,1)./result_p(2:end,1));
result_order_p = table(result_order_p(:,1),result_order_p(:,2),result_order_p(:,3),result_order_p(:,4),'VariableNames',{'Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
draw(f_1,f_2,result_p,['(Linear 2D FEM)',"h","h^2"]);
result_p = table(sym(result_p(:,1)),sym(result_p(:,2)),sym(result_p(:,3)),result_p(:,4),result_p(:,5),result_p(:,6),result_p(:,7),...
               'VariableNames',{'hx','hy','ht','Max Absolute Error','L^infty Error','L^2 Error','H^1(semi) Error'});
fprintf("p: \n");
disp(result_p);
fprintf("Convergence Order: \n");
disp(result_order_p);
%% 绘图函数
function draw(f_1,f_2,data,txt)
x = data(:,1);
% 绘制 O(*) 的参考线
figure;
loglog(x,f_1(x),'-o','LineWidth',2,'Color','r');
hold on;
loglog(x,f_2(x),'-o','LineWidth',2,'Color','g');
% 绘制实际的误差图
loglog(x,data(:,4), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,5), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,6), '--','LineWidth',2,'Marker',"pentagram");
loglog(x,data(:,7), '--','LineWidth',2,'Marker',"pentagram");
grid on;
grid minor;
xlabel('log10(h)');
ylabel('log10(error)');
title(string(txt(1)) + ' Error vs Mesh Width');
legend("$O(" + txt(2) +  ")$ reference","$O(" + txt(3) +  ")$ reference",...
        'max abs error for all finite nodes','$L^{\infty}$ error','$L^{2}$ error','$H^1(semi)$ error',...
        'Interpreter','latex','Location','best');
end